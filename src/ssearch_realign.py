"""
Module for rigorous realignment using SSEARCH (Smith-Waterman).
Implements pipeline step 3: Rigorous Realignment with statistical permutation.
"""

import subprocess
import logging
import tempfile
import multiprocessing
from pathlib import Path
from typing import List, Optional, Dict, Any, Tuple
from dataclasses import dataclass
from functools import partial

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.config import SSearchParameters, DEFAULT_SSEARCH_PARAMS, SSEARCH_EVALUE_THRESHOLD
from src.blast_filter import BlastHit


logger = logging.getLogger(__name__)


class SSearchError(Exception):
    """Exception for SSEARCH execution errors."""
    pass


@dataclass
class SSearchAlignment:
    """Represents an SSEARCH alignment."""
    
    query_id: str
    subject_id: str
    identity: float
    alignment_length: int
    mismatches: int
    gap_opens: int
    query_start: int
    query_end: int
    subject_start: int
    subject_end: int
    evalue: float
    bit_score: float
    raw_score: int
    query_seq: str
    subject_seq: str
    
    @classmethod
    def from_m10_line(cls, line: str) -> Optional['SSearchAlignment']:
        """
        Create SSearchAlignment object from format 8C line.
        
        Format 8C (12 fields, same as BLAST): qseqid sseqid pident length 
                    mismatch gapopen qstart qend sstart send evalue bitscore
        """
        fields = line.strip().split('\t')
        if len(fields) < 12:
            return None
        
        try:
            return cls(
                query_id=fields[0],
                subject_id=fields[1],
                identity=float(fields[2]),
                alignment_length=int(fields[3]),
                mismatches=int(fields[4]),
                gap_opens=int(fields[5]),
                query_start=int(fields[6]),
                query_end=int(fields[7]),
                subject_start=int(fields[8]),
                subject_end=int(fields[9]),
                evalue=float(fields[10]),
                bit_score=float(fields[11]),
                raw_score=0,  # Será preenchido se disponível
                query_seq="",
                subject_seq=""
            )
        except (ValueError, IndexError) as e:
            logger.warning(f"Error parsing SSEARCH line: {line[:50]}... - {e}")
            return None


def create_fasta_from_sequence(
    sequence: str,
    seq_id: str,
    output_path: Path,
    description: str = ""
) -> Path:
    """
    Create temporary FASTA file from a sequence.
    
    Args:
        sequence: Sequence to write
        seq_id: Sequence identifier
        output_path: Path to output file
        description: Optional sequence description
    
    Returns:
        Path to created FASTA file
    """
    record = SeqRecord(
        Seq(sequence),
        id=seq_id,
        description=description
    )
    
    with open(output_path, 'w') as f:
        SeqIO.write(record, f, "fasta")
    
    return output_path


def run_ssearch36(
    query_fasta: Path,
    target_fasta: Path,
    output_file: Path,
    ssearch_path: str,
    ssearch_params: SSearchParameters = DEFAULT_SSEARCH_PARAMS
) -> Path:
    """
    Execute SSEARCH36 (Smith-Waterman) with statistical permutation.
    
    Args:
        query_fasta: FASTA file with query sequence (functional protein)
        target_fasta: FASTA file with target sequence (genomic)
        output_file: Output file for results
        ssearch_path: Path to ssearch36 executable
        ssearch_params: SSEARCH parameters
    
    Returns:
        Path to output file
    
    Raises:
        SSearchError: If execution fails
    """
    if not query_fasta.exists():
        raise SSearchError(f"Arquivo query não encontrado: {query_fasta}")
    
    if not target_fasta.exists():
        raise SSearchError(f"Arquivo target não encontrado: {target_fasta}")
    
    logger.debug(f"Running SSEARCH36: {query_fasta.name} vs {target_fasta.name}")
    
    # Build command with full ssearch36 path
    command = [ssearch_path]
    command.extend(ssearch_params.to_command_args())
    command.extend([str(query_fasta), str(target_fasta)])
    
    logger.debug(f"SSEARCH command: {' '.join(command)}")
    
    try:
        with open(output_file, 'w') as out:
            result = subprocess.run(
                command,
                stdout=out,
                stderr=subprocess.PIPE,  # Silence SSEARCH messages
                text=True,
                check=True
            )
        
        logger.debug(f"SSEARCH complete. Results in: {output_file}")
        return output_file
        
    except subprocess.CalledProcessError as e:
        error_msg = f"Erro ao executar ssearch36: {e.stderr}"
        logger.error(error_msg)
        raise SSearchError(error_msg) from e
    except FileNotFoundError:
        error_msg = "ssearch36 não encontrado. Certifique-se de que o FASTA36 está instalado."
        logger.error(error_msg)
        raise SSearchError(error_msg)


def parse_ssearch_output(ssearch_output_path: Path) -> List[SSearchAlignment]:
    """
    Parse SSEARCH output file in m10 (tabular) format.
    
    Args:
        ssearch_output_path: Path to SSEARCH output file
    
    Returns:
        List of SSearchAlignment objects
    """
    logger.debug(f"Fazendo parsing do arquivo SSEARCH: {ssearch_output_path}")
    
    if not ssearch_output_path.exists():
        raise FileNotFoundError(f"Arquivo SSEARCH não encontrado: {ssearch_output_path}")
    
    alignments = []
    
    with open(ssearch_output_path, 'r') as f:
        for line in f:
            # Ignore comments
            if line.startswith('#'):
                continue
            
            alignment = SSearchAlignment.from_m10_line(line)
            if alignment:
                alignments.append(alignment)
    
    logger.debug(f"Total SSEARCH alignments read: {len(alignments)}")
    return alignments


def _ssearch_worker(
    args: Tuple[BlastHit, Path, Path, Path, str, SSearchParameters]
) -> Tuple[str, Optional[SSearchAlignment]]:
    """
    Isolated worker for parallel SSEARCH execution.
    
    Args:
        args: Tuple with (blast_hit, query_fasta_path, genome_fasta_path, 
              work_dir, ssearch_path, ssearch_params)
    
    Returns:
        Tuple (hit_key, alignment) where alignment can be None if it fails
    """
    blast_hit, query_fasta_path, genome_fasta_path, work_dir, ssearch_path, ssearch_params = args
    
    hit_key = f"{blast_hit.qseqid}_{blast_hit.sseqid}_{blast_hit.sstart}_{blast_hit.send}"
    
    alignment = realign_blast_hit_with_ssearch(
        blast_hit,
        query_fasta_path,
        genome_fasta_path,
        work_dir,
        ssearch_path,
        ssearch_params
    )
    
    return (hit_key, alignment)


def realign_blast_hit_with_ssearch(
    blast_hit: BlastHit,
    query_fasta_path: Path,
    genome_fasta_path: Path,
    work_dir: Path,
    ssearch_path: str,
    ssearch_params: SSearchParameters = DEFAULT_SSEARCH_PARAMS
) -> Optional[SSearchAlignment]:
    """
    Realign a BLAST hit using SSEARCH with statistical permutation.
    
    Args:
        blast_hit: BLAST hit to be realigned
        query_fasta_path: FASTA file with original query sequences
        genome_fasta_path: Genome FASTA file
        work_dir: Working directory for temporary files
        ssearch_path: Path to ssearch36 executable
        ssearch_params: SSEARCH parameters
    
    Returns:
        SSearchAlignment object or None if it fails
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    
    # Create temporary FASTA file with query sequence
    query_temp = work_dir / f"query_{blast_hit.qseqid}_{blast_hit.sstart}.fasta"
    create_fasta_from_sequence(
        blast_hit.qseq.replace('-', ''),  # Remove gaps
        blast_hit.qseqid,
        query_temp
    )
    
    # Create temporary FASTA file with subject sequence (genomic)
    target_temp = work_dir / f"target_{blast_hit.sseqid}_{blast_hit.sstart}.fasta"
    create_fasta_from_sequence(
        blast_hit.sseq.replace('-', ''),  # Remove gaps
        f"{blast_hit.sseqid}:{blast_hit.sstart}-{blast_hit.send}",
        target_temp
    )
    
    # Output file
    output_file = work_dir / f"ssearch_{blast_hit.qseqid}_{blast_hit.sstart}.out"
    
    try:
        # Execute SSEARCH
        run_ssearch36(query_temp, target_temp, output_file, ssearch_path, ssearch_params)
        
        # Parse result
        alignments = parse_ssearch_output(output_file)
        
        if alignments:
            return alignments[0]  # Return best alignment
        else:
            logger.debug(f"No SSEARCH alignment found for hit: {blast_hit.qseqid}")
            return None
            
    except SSearchError as e:
        logger.error(f"Error realigning hit with SSEARCH: {e}")
        return None
    finally:
        # Cleanup temporary files (optional)
        # query_temp.unlink(missing_ok=True)
        # target_temp.unlink(missing_ok=True)
        pass


def realign_filtered_hits(
    filtered_hits: List[BlastHit],
    query_fasta_path: Path,
    genome_fasta_path: Path,
    output_dir: Path,
    ssearch_path: str,
    ssearch_params: SSearchParameters = DEFAULT_SSEARCH_PARAMS
) -> Dict[str, SSearchAlignment]:
    """
    Realign all filtered hits using SSEARCH.
    
    Args:
        filtered_hits: List of filtered BLAST hits
        query_fasta_path: FASTA file with query sequences
        genome_fasta_path: Genome FASTA file
        output_dir: Output directory for files
        ssearch_path: Path to ssearch36 executable
        ssearch_params: SSEARCH parameters
    
    Returns:
        Dictionary mapping hit identifier to its SSEARCH alignment
    """
    logger.debug(f"Realigning {len(filtered_hits)} hits with SSEARCH")
    
    ssearch_work_dir = output_dir / "ssearch_work"
    ssearch_work_dir.mkdir(parents=True, exist_ok=True)
    
    realignments = {}
    
    for i, hit in enumerate(filtered_hits, 1):
        logger.debug(f"Realigning hit {i}/{len(filtered_hits)}: {hit.qseqid}")
        
        alignment = realign_blast_hit_with_ssearch(
            hit,
            query_fasta_path,
            genome_fasta_path,
            ssearch_work_dir,
            ssearch_path,
            ssearch_params
        )
        
        if alignment:
            # Filter 2: Apply E-value threshold (1e-6)
            if alignment.evalue <= SSEARCH_EVALUE_THRESHOLD:
                hit_key = f"{hit.qseqid}_{hit.sseqid}_{hit.sstart}_{hit.send}"
                realignments[hit_key] = alignment
            else:
                logger.debug(f"SSEARCH alignment filtered by E-value: {alignment.evalue:.2e} > {SSEARCH_EVALUE_THRESHOLD:.2e}")
    
    logger.debug(f"Realignment complete. Total: {len(realignments)} alignments (after E-value filter)")
    return realignments


def realign_filtered_hits_parallel(
    filtered_hits: List[BlastHit],
    query_fasta_path: Path,
    genome_fasta_path: Path,
    output_dir: Path,
    ssearch_path: str,
    num_threads: Optional[int] = None,
    ssearch_params: SSearchParameters = DEFAULT_SSEARCH_PARAMS
) -> Dict[str, SSearchAlignment]:
    """
    Realign all filtered hits using SSEARCH in parallel.
    
    Args:
        filtered_hits: List of filtered BLAST hits
        query_fasta_path: FASTA file with query sequences
        genome_fasta_path: Genome FASTA file
        output_dir: Output directory for files
        ssearch_path: Path to ssearch36 executable
        num_threads: Number of threads (None = CPU_COUNT-1)
        ssearch_params: SSEARCH parameters
    
    Returns:
        Dictionary mapping hit identifier to its SSEARCH alignment
    """

    output_dir.mkdir(parents=True, exist_ok=True)

    # Detect number of CPUs if not specified
    if num_threads is None:
        num_threads = max(1, multiprocessing.cpu_count() - 1)
    
    logger.debug(f"Realigning {len(filtered_hits)} hits with SSEARCH in parallel")
    logger.debug(f"Using {num_threads} threads")
    
    ssearch_work_dir = output_dir / "ssearch_work"
    ssearch_work_dir.mkdir(parents=True, exist_ok=True)
    
    # Prepare arguments for each worker
    worker_args = [
        (hit, query_fasta_path, genome_fasta_path, ssearch_work_dir, ssearch_path, ssearch_params)
        for hit in filtered_hits
    ]
    
    realignments = {}
    
    # Execute in parallel
    with multiprocessing.Pool(processes=num_threads) as pool:
        results = pool.map(_ssearch_worker, worker_args)
    
    # Process results
    for hit_key, alignment in results:
        if alignment:
            # Filter 2: Apply E-value threshold (1e-6)
            if alignment.evalue <= SSEARCH_EVALUE_THRESHOLD:
                realignments[hit_key] = alignment
            else:
                logger.debug(f"SSEARCH alignment filtered by E-value: {alignment.evalue:.2e} > {SSEARCH_EVALUE_THRESHOLD:.2e}")

    realignments_output = output_dir / "ssearch_realignments.tsv"
    with open(realignments_output, 'w') as f:
        f.write("hit_key\tquery_id\tsubject_id\tidentity\talignment_length\tmismatches\tgap_opens\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbit_score\traw_score\tquery_seq\tsubject_seq\n")
        for key, aln in realignments.items():
            f.write(f"{key}\t{aln.query_id}\t{aln.subject_id}\t{aln.identity:.2f}\t"
                    f"{aln.alignment_length}\t{aln.mismatches}\t{aln.gap_opens}\t"
                    f"{aln.query_start}\t{aln.query_end}\t{aln.subject_start}\t"
                    f"{aln.subject_end}\t{aln.evalue:.2e}\t{aln.bit_score:.2f}\t"
                    f"{aln.raw_score}\t{aln.query_seq}\t{aln.subject_seq}\n")
        
        # for key, aln in realignments.items():
        #     f.write(f"{key}\t{aln.query_id}\t{aln.subject_id}\t"
        #            f"{aln.identity:.2f}\t{aln.evalue:.2e}\t{aln.bit_score:.2f}\n")
    
    logger.debug(f"Parallel realignment complete. Total: {len(realignments)} alignments (after E-value filter)")
    return realignments


import pandas as pd
import sys
def load_blast_hits_from_tsv(tsv_path: Path) -> List[BlastHit]:
    """
    Args:
        tsv_path: Path to tsv containing hits
    
    Returns:
        List containing BlastHit objects
    """
    # Lê o TSV
    df = pd.read_csv(tsv_path, sep='\t')
    
    # Converte o DataFrame de volta para uma lista de objetos BlastHit
    # O to_dict('records') cria uma lista de dicionários onde as chaves são as colunas
    hits = [BlastHit(**row) for row in df.to_dict('records')]
    return hits

if __name__ == "__main__":
    tsv_file = sys.argv[1]
    hits = load_blast_hits_from_tsv(tsv_file)
    realign_filtered_hits_parallel(hits, Path(sys.argv[2]), Path(sys.argv[3]), Path(sys.argv[4]), sys.argv[5])