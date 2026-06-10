"""
Module for BLAST database preparation and search execution.
Implements pipeline step 1: Preparation and Homology Search.
"""

import subprocess
import logging
from pathlib import Path
from typing import Optional

from src.config import BlastParameters, DEFAULT_BLAST_PARAMS


logger = logging.getLogger(__name__)


class BlastDatabaseError(Exception):
    """Exception for BLAST database creation errors."""
    pass


class BlastSearchError(Exception):
    """Exception for BLAST search execution errors."""
    pass


def create_blast_database(
    genome_fasta_path: Path,
    output_db_path: Optional[Path] = None,
    dbtype: str = "nucl",
    makeblastdb_path: str = "makeblastdb"
) -> Path:
    """
    Create a BLAST database from a FASTA file.
    
    Args:
        genome_fasta_path: Path to genome FASTA file
        output_db_path: Path to output database (optional)
        dbtype: Database type ('nucl' for nucleotides, 'prot' for proteins)
    
    Returns:
        Path to created database
    
    Raises:
        BlastDatabaseError: If database creation fails
    """
    if not genome_fasta_path.exists():
        raise BlastDatabaseError(f"Arquivo de genoma não encontrado: {genome_fasta_path}")
    
    # Use same path as input file if not specified
    if output_db_path is None:
        output_db_path = genome_fasta_path
    
    logger.debug(f"Creating BLAST database...")
    
    command = [
        makeblastdb_path,
        "-in", str(genome_fasta_path),
        "-dbtype", dbtype,
        "-out", str(output_db_path)
    ]
    
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True
        )
        logger.debug("BLAST database created")
        return output_db_path
        
    except subprocess.CalledProcessError as e:
        error_msg = f"Erro ao criar banco de dados BLAST: {e.stderr}"
        logger.error(error_msg)
        raise BlastDatabaseError(error_msg) from e
    except FileNotFoundError:
        error_msg = "makeblastdb não encontrado. Certifique-se de que o BLAST+ está instalado."
        logger.error(error_msg)
        raise BlastDatabaseError(error_msg)


def run_tblastn_search(
    query_protein_fasta: Path,
    database_path: Path,
    output_file: Path,
    blast_params: BlastParameters = DEFAULT_BLAST_PARAMS,
    tblastn_path: str = "tblastn"
) -> Path:
    """
    Execute tblastn search with strict parameters.
    
    Search protein sequences (query) against nucleotide database,
    translating database in all 6 reading frames.
    
    Args:
        query_protein_fasta: FASTA file with query protein sequences (lysozymes)
        database_path: Path to nucleotide BLAST database
        output_file: Output file for results
        blast_params: BLAST parameters (uses DEFAULT_BLAST_PARAMS if not specified)
    
    Returns:
        Path to output file with results
    
    Raises:
        BlastSearchError: If search execution fails
    """
    if not query_protein_fasta.exists():
        raise BlastSearchError(f"Arquivo query não encontrado: {query_protein_fasta}")
    
    # Create output directory if it doesn't exist
    output_file.parent.mkdir(parents=True, exist_ok=True)
    
    logger.debug(f"Running tblastn with {blast_params.num_threads} threads...")
    
    # Build command with parameters
    command = [
        tblastn_path,
        "-query", str(query_protein_fasta),
        "-db", str(database_path),
        "-out", str(output_file)
    ]
    
    # Adiciona os parâmetros customizados
    command.extend(blast_params.to_command_args())
    
    logger.debug(f"BLAST command: {' '.join(command)}")
    
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            check=True
        )
        
        logger.debug("BLAST search complete")
        return output_file
        
    except subprocess.CalledProcessError as e:
        # Check if it's a segmentation fault (SIGSEGV)
        if e.returncode == -11 or 'SIGSEGV' in str(e):
            logger.warning("BLAST crashed with SIGSEGV. Attempting retry with reduced threads...")
            
            # Retry with single thread as workaround for BLAST segfault bug
            retry_command = command.copy()
            for i, arg in enumerate(retry_command):
                if arg == '-num_threads':
                    retry_command[i+1] = '1'
                    break
            
            try:
                logger.debug(f"Retry command: {' '.join(retry_command)}")
                result = subprocess.run(
                    retry_command,
                    capture_output=True,
                    text=True,
                    check=True,
                    timeout=600  # 10 minute timeout
                )
                logger.info("BLAST succeeded with single thread (workaround for segfault)")
                return output_file
                
            except (subprocess.CalledProcessError, subprocess.TimeoutExpired) as retry_error:
                error_msg = (
                    f"BLAST failed with SIGSEGV (segmentation fault).\n"
                    f"This usually indicates:\n"
                    f"  1. Corrupted BLAST+ installation\n"
                    f"  2. Incompatible BLAST+ version\n"
                    f"  3. Memory corruption or hardware issues\n\n"
                    f"Troubleshooting steps:\n"
                    f"  - Reinstall BLAST+: sudo apt remove ncbi-blast+ && sudo apt install ncbi-blast+\n"
                    f"  - Or download latest from NCBI: https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/\n"
                    f"  - Check system resources: free -h\n"
                    f"  - Reduce threads: --num-threads 1\n\n"
                    f"Original error: {e.stderr if e.stderr else 'SIGSEGV'}\n"
                    f"Retry error: {retry_error}"
                )
                logger.error(error_msg)
                raise BlastSearchError(error_msg) from e
        else:
            error_msg = f"Erro ao executar tblastn: {e.stderr}"
            logger.error(error_msg)
            raise BlastSearchError(error_msg) from e
            
    except FileNotFoundError:
        error_msg = "tblastn não encontrado. Certifique-se de que o BLAST+ está instalado."
        logger.error(error_msg)
        raise BlastSearchError(error_msg)
    except subprocess.TimeoutExpired:
        error_msg = "BLAST search timeout (>10 minutes). Check system resources or reduce dataset size."
        logger.error(error_msg)
        raise BlastSearchError(error_msg)


def run_blast_pipeline_step(
    genome_fasta: Path,
    lysozyme_fasta: Path,
    output_dir: Path,
    blast_params: BlastParameters = DEFAULT_BLAST_PARAMS,
    makeblastdb_path: str = "makeblastdb",
    tblastn_path: str = "tblastn"
) -> Path:
    """
    Execute complete BLAST preparation and search step.
    
    Args:
        genome_fasta: E. coli genome FASTA file
        lysozyme_fasta: Reference lysozyme FASTA file
        output_dir: Output directory for files
        blast_params: BLAST parameters
    
    Returns:
        Path to file with BLAST results
    """
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create BLAST database
    db_path = output_dir / "genome_db"
    
    create_blast_database(genome_fasta, db_path, makeblastdb_path=makeblastdb_path)

    # Execute tblastn search
    blast_output = output_dir / "blast_results.tsv"
    
    run_tblastn_search(lysozyme_fasta, db_path, blast_output, blast_params, tblastn_path)
    return blast_output

import sys
if __name__ == '__main__':
    run_blast_pipeline_step(Path(sys.argv[1]), Path(sys.argv[2]), Path(sys.argv[3]))