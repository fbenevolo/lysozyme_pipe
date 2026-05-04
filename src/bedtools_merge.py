"""
Module for merging adjacent segments using BEDTools.
Implements pipeline step 4: Segment Merging.
"""

import subprocess
import logging
import pandas as pd
from pathlib import Path
from typing import List, Dict, Any
from dataclasses import dataclass

from src.config import BedToolsParameters, DEFAULT_BEDTOOLS_PARAMS
from src.blast_filter import BlastHit


logger = logging.getLogger(__name__)


class BedToolsError(Exception):
    """Exception for BEDTools execution errors."""
    pass


@dataclass
class GenomicRegion:
    """Represents a merged genomic region."""
    
    genome_id: str          # Identificador do genoma de origem
    chromosome: str         # Sequência/cromossomo
    start: int              # Posição inicial
    end: int                # Posição final
    strand: str             # Fita (+/-)
    num_hsps: int           # Número de HSPs fundidos
    mean_score: float       # Score médio dos HSPs
    min_score: float        # Score mínimo
    max_score: float        # Score máximo
    query_ids: List[str]    # IDs das queries que mapearam nesta região
    
    @property
    def length(self) -> int:
        """Return region length."""
        return self.end - self.start + 1
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert region to dictionary."""
        return {
            'genome_id': self.genome_id,
            'chromosome': self.chromosome,
            'start': self.start,
            'end': self.end,
            'strand': self.strand,
            'length': self.length,
            'num_hsps': self.num_hsps,
            'mean_score': self.mean_score,
            'min_score': self.min_score,
            'max_score': self.max_score,
            'query_ids': ','.join(self.query_ids) if isinstance(self.query_ids, list) else self.query_ids
        }


def blast_hits_to_bed(
    hits: List[BlastHit],
    output_bed_path: Path
) -> Path:
    """
    Convert BLAST hits to BED format.
    
    BED format: chromosome start end name score strand
    
    Args:
        hits: List of BLAST hits
        output_bed_path: Path to output BED file
    
    Returns:
        Path to created BED file
    """
    logger.debug(f"Converting hits to BED format...")
    
    bed_entries = []
    
    for hit in hits:
        # Ensure start < end (BED uses 0-based coordinates)
        start = min(hit.sstart, hit.send) - 1  # Convert to 0-based
        end = max(hit.sstart, hit.send)
        
        # Determine strand
        strand = "+" if hit.sstart < hit.send else "-"
        
        bed_entry = {
            'chromosome': hit.sseqid,
            'start': start,
            'end': end,
            'name': hit.qseqid,
            'score': hit.score,
            'strand': strand
        }
        bed_entries.append(bed_entry)
    
    # Create DataFrame and sort
    df = pd.DataFrame(bed_entries)
    df = df.sort_values(['chromosome', 'start', 'end'])
    
    # Save in BED format
    output_bed_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(
        output_bed_path,
        sep='\t',
        header=False,
        index=False
    )
    
    logger.debug(f"BED file created with {len(df)} entries")
    return output_bed_path


def run_bedtools_merge(
    input_bed_path: Path,
    output_bed_path: Path,
    bedtools_params: BedToolsParameters = DEFAULT_BEDTOOLS_PARAMS,
    bedtools_path: str = "bedtools"
) -> Path:
    """
    Execute bedtools merge to merge overlapping or adjacent regions.
    
    Args:
        input_bed_path: Input BED file
        output_bed_path: Output BED file
        bedtools_params: BEDTools parameters
    
    Returns:
        Path to merged BED file
    
    Raises:
        BedToolsError: If execution fails
    """
    logger.debug(f"Running bedtools merge...")
    
    if not input_bed_path.exists():
        raise BedToolsError(f"Arquivo BED não encontrado: {input_bed_path}")
    
    # Build command
    command = [bedtools_path, "merge", "-i", str(input_bed_path)]
    command.extend(bedtools_params.to_command_args())
    
    logger.debug(f"BEDTools command: {' '.join(command)}")
    
    try:
        with open(output_bed_path, 'w') as out:
            result = subprocess.run(
                command,
                stdout=out,
                stderr=subprocess.PIPE,
                text=True,
                check=True
            )
        
        logger.debug(f"bedtools merge complete")
        return output_bed_path
        
    except subprocess.CalledProcessError as e:
        error_msg = f"Erro ao executar bedtools merge: {e.stderr}"
        logger.error(error_msg)
        raise BedToolsError(error_msg) from e
    except FileNotFoundError:
        error_msg = "bedtools não encontrado. Certifique-se de que o BEDTools está instalado."
        logger.error(error_msg)
        raise BedToolsError(error_msg)


def parse_merged_bed(merged_bed_path: Path, genome_id: str = "unknown") -> List[GenomicRegion]:
    """
    Parse merged BED file generated by bedtools merge.
    
    bedtools merge output file with -c 2,4,4,4,5 -o count,mean,min,max,collapse has:
    - Columns 1-3: chromosome, start, end
    - Column 4: count of names (HSPs)
    - Column 5: mean score
    - Column 6: min score
    - Column 7: max score
    - Column 8: collapsed strands
    
    Args:
        merged_bed_path: Path to merged BED file
        genome_id: Genome identifier
    
    Returns:
        List of GenomicRegion objects
    """
    logger.debug(f"Parsing merged BED file: {merged_bed_path}")
    
    if not merged_bed_path.exists():
        raise FileNotFoundError(f"Arquivo BED fundido não encontrado: {merged_bed_path}")
    
    # Read file
    df = pd.read_csv(
        merged_bed_path,
        sep='\t',
        header=None,
        names=['chromosome', 'start', 'end', 'num_hsps', 'mean_score', 
               'min_score', 'max_score', 'strands']
    )
    
    regions = []
    
    for _, row in df.iterrows():
        # Determine predominant strand
        strands = str(row['strands']).split(',')
        strand = max(set(strands), key=strands.count) if strands else '+'
        
        # Treat '.' values as 0.0 (when bedtools cannot calculate)
        try:
            mean_score = float(row['mean_score']) if row['mean_score'] != '.' else 0.0
            min_score = float(row['min_score']) if row['min_score'] != '.' else 0.0
            max_score = float(row['max_score']) if row['max_score'] != '.' else 0.0
        except (ValueError, TypeError):
            mean_score = min_score = max_score = 0.0
        
        region = GenomicRegion(
            genome_id=genome_id,
            chromosome=row['chromosome'],
            start=int(row['start']),
            end=int(row['end']),
            strand=strand,
            num_hsps=int(row['num_hsps']),
            mean_score=mean_score,
            min_score=min_score,
            max_score=max_score,
            query_ids=[]  # Will be populated later
        )
        regions.append(region)
    
    logger.debug(f"Total merged regions: {len(regions)}")
    return regions


def merge_blast_hits(
    hits: List[BlastHit],
    output_dir: Path,
    bedtools_params: BedToolsParameters = DEFAULT_BEDTOOLS_PARAMS,
    bedtools_path: str = "bedtools",
    genome_id: str = "unknown"
) -> List[GenomicRegion]:
    """
    Execute complete BLAST hit merging step.
    
    Args:
        hits: List of filtered BLAST hits
        output_dir: Output directory for intermediate and output files
        bedtools_params: BEDTools parameters
        bedtools_path: Path to bedtools executable
        genome_id: Genome identifier
    
    Returns:
        List of merged genomic regions
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    
    # Convert hits to BED
    bed_path = output_dir / "blast_hits.bed"
    blast_hits_to_bed(hits, bed_path)
    
    # Execute bedtools merge
    merged_bed_path = output_dir / "merged_regions.bed"
    run_bedtools_merge(bed_path, merged_bed_path, bedtools_params, bedtools_path)
    
    # Parse merged regions
    regions = parse_merged_bed(merged_bed_path, genome_id=genome_id)
    
    return regions
    

import sys
from src.ssearch_realign import load_hits_from_tsv
if __name__ == "__main__":
    blast_hits = load_hits_from_tsv(sys.argv[1])

    logger.info(sys.argv[3])

    merge_blast_hits(blast_hits, Path(sys.argv[2]), genome_id=sys.argv[3])