"""
Module for parsing and filtering BLAST results.
Implements pipeline step 2: Initial Filtering (Quality Control).
"""

import pandas as pd
import logging
from pathlib import Path
from typing import List, Dict, Any
from dataclasses import dataclass

from src.config import MIN_IDENTITY_THRESHOLD, MIN_BLOSUM62_SCORE


logger = logging.getLogger(__name__)


@dataclass
class BlastHit:
    """Represents an individual BLAST hit."""
    
    genome_id: str   # Identificador do genoma de origem
    qseqid: str      # Query sequence ID
    qlen: int        # Query length
    sseqid: str      # Subject sequence ID
    slen: int        # Subject length
    qstart: int      # Query start position
    qend: int        # Query end position
    sstart: int      # Subject start position
    send: int        # Subject end position
    qseq: str        # Query aligned sequence
    sseq: str        # Subject aligned sequence
    evalue: float    # E-value
    bitscore: float  # Bit score
    score: int       # Raw score (BLOSUM62)
    length: int      # Alignment length
    pident: float    # Percent identity
    nident: int      # Number of identical matches
    mismatch: int    # Number of mismatches
    positive: int    # Number of positive matches
    gapopen: int     # Number of gap openings
    gaps: int        # Total gaps
    ppos: float      # Percent positive matches
    sframe: int      # Subject frame
    sstrand: str     # Subject strand
    qcovs: float     # Query coverage per subject
    qcovhsp: float   # Query coverage per HSP
    
    def passes_quality_filter(
        self,
        min_identity: float = MIN_IDENTITY_THRESHOLD,
        min_score: int = MIN_BLOSUM62_SCORE
    ) -> bool:
        """
        Check if hit passes quality criteria.
        
        Discard rule: Remove alignments that have ANY of these conditions:
        - %identity <= min_identity (default: 20%)
        - BLOSUM62 score < min_score (default: 120)
        
        Args:
            min_identity: Minimum percent identity threshold
            min_score: Minimum BLOSUM62 score
        
        Returns:
            True if hit passes filter (should be kept), False otherwise
        """
        # Discard if ANY condition is true (OR logic)
        fails_identity = self.pident <= min_identity
        fails_score = self.score < min_score
        
        if fails_identity or fails_score:
            logger.debug(
                f"Hit discarded - qseqid: {self.qseqid}, sseqid: {self.sseqid}, "
                f"pident: {self.pident}%, score: {self.score}"
            )
            return False
        
        return True
    
    def to_dict(self) -> Dict[str, Any]:
        """Convert hit to dictionary."""
        return {
            'genome_id': self.genome_id,
            'qseqid': self.qseqid,
            'qlen': self.qlen,
            'sseqid': self.sseqid,
            'slen': self.slen,
            'qstart': self.qstart,
            'qend': self.qend,
            'sstart': self.sstart,
            'send': self.send,
            'qseq': self.qseq,
            'sseq': self.sseq,
            'evalue': self.evalue,
            'bitscore': self.bitscore,
            'score': self.score,
            'length': self.length,
            'pident': self.pident,
            'nident': self.nident,
            'mismatch': self.mismatch,
            'positive': self.positive,
            'gapopen': self.gapopen,
            'gaps': self.gaps,
            'ppos': self.ppos,
            'sframe': self.sframe,
            'sstrand': self.sstrand,
            'qcovs': self.qcovs,
            'qcovhsp': self.qcovhsp
        }
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'BlastHit':
        return cls(**d)


def parse_blast_output(blast_output_path: Path, genome_id: str = "unknown") -> List[BlastHit]:
    """
    Parse BLAST output file in tabular format.
    
    Args:
        blast_output_path: Path to BLAST output file
        genome_id: Genome identifier
    
    Returns:
        List of BlastHit objects
    """
    logger.debug(f"Parsing BLAST file: {blast_output_path}")
    
    if not blast_output_path.exists():
        raise FileNotFoundError(f"Arquivo BLAST não encontrado: {blast_output_path}")
    
    # Define column names as specified in config.py
    column_names = [
        "qseqid", "qlen", "sseqid", "slen", "qstart", "qend", "sstart", "send",
        "qseq", "sseq", "evalue", "bitscore", "score", "length", "pident",
        "nident", "mismatch", "positive", "gapopen", "gaps", "ppos", "sframe",
        "sstrand", "qcovs", "qcovhsp"
    ]
    
    # Read tabular file
    df = pd.read_csv(
        blast_output_path,
        sep="\t",
        names=column_names,
        comment="#"
    )
    
    logger.debug(f"Total BLAST hits read: {len(df)}")
    
    # Convert each row to BlastHit object
    hits = []
    for _, row in df.iterrows():
        hit = BlastHit(
            genome_id=genome_id,
            qseqid=row['qseqid'],
            qlen=int(row['qlen']),
            sseqid=row['sseqid'],
            slen=int(row['slen']),
            qstart=int(row['qstart']),
            qend=int(row['qend']),
            sstart=int(row['sstart']),
            send=int(row['send']),
            qseq=row['qseq'],
            sseq=row['sseq'],
            evalue=float(row['evalue']),
            bitscore=float(row['bitscore']),
            score=int(row['score']),
            length=int(row['length']),
            pident=float(row['pident']),
            nident=int(row['nident']),
            mismatch=int(row['mismatch']),
            positive=int(row['positive']),
            gapopen=int(row['gapopen']),
            gaps=int(row['gaps']),
            ppos=float(row['ppos']),
            sframe=int(row['sframe']),
            sstrand=row['sstrand'],
            qcovs=float(row['qcovs']),
            qcovhsp=float(row['qcovhsp'])
        )
        hits.append(hit)
    
    return hits


def filter_blast_hits(
    hits: List[BlastHit],
    min_identity: float = MIN_IDENTITY_THRESHOLD,
    min_score: int = MIN_BLOSUM62_SCORE
) -> List[BlastHit]:
    """
    Filter BLAST hits based on quality criteria.
    
    Remove alignments that have BOTH conditions:
    - %identity <= min_identity
    - BLOSUM62 score < min_score
    
    Args:
        hits: List of BLAST hits
        min_identity: Minimum percent identity threshold
        min_score: Minimum BLOSUM62 score
    
    Returns:
        List of filtered hits
    """
    logger.info(f"Filtering {len(hits)} BLAST hits...")
    
    filtered_hits = [
        hit for hit in hits
        if hit.passes_quality_filter(min_identity, min_score)
    ]
    
    removed_count = len(hits) - len(filtered_hits)
    logger.info(f"Kept: {len(filtered_hits)}, Removed: {removed_count}")
    
    return filtered_hits


def save_filtered_hits(hits: List[BlastHit], output_path: Path) -> None:
    """
    Save filtered hits to TSV file.
    
    Args:
        hits: List of BLAST hits
        output_path: Path to output file
    """
    logger.debug(f"Saving {len(hits)} filtered hits to: {output_path}")
    
    # Convert to DataFrame
    data = [hit.to_dict() for hit in hits]
    df = pd.DataFrame(data)
    
    # Save as TSV
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep="\t", index=False)
    
    logger.debug("Filtered hits saved successfully")


def run_filtering_step(
    blast_output_path: Path,
    filtered_output_path: Path,
    min_identity: float = MIN_IDENTITY_THRESHOLD,
    min_score: int = MIN_BLOSUM62_SCORE,
    genome_id: str = "unknown"
) -> List[BlastHit]:
    """
    Execute complete BLAST hit filtering step.
    
    Args:
        blast_output_path: Path to BLAST output file
        filtered_output_path: Path to save filtered hits
        min_identity: Minimum percent identity threshold
        min_score: Minimum BLOSUM62 score
        genome_id: Genome identifier
    
    Returns:
        List of filtered hits
    """
    # Parse BLAST file with genome_id
    hits = parse_blast_output(blast_output_path, genome_id=genome_id)

    # Filtering
    filtered_hits = filter_blast_hits(hits, min_identity, min_score)
    
    # Save results
    save_filtered_hits(filtered_hits, filtered_output_path)
    
    return filtered_hits


import sys
if __name__ == '__main__':
    run_filtering_step(Path(sys.argv[1]), Path(sys.argv[2]), genome_id=sys.argv[3])