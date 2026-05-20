"""
Module for score density calculation and best hit selection.
Implements pipeline step 5: Score Density Calculation and Selection.
"""

import pandas as pd
import logging
from pathlib import Path
from typing import List, Dict, Tuple, Optional
from dataclasses import dataclass
from collections import defaultdict

from src.blast_filter import BlastHit
from src.bedtools_merge import GenomicRegion


logger = logging.getLogger(__name__)


@dataclass
class ProteinHitGroup:
    """Group of HSPs for a specific protein in a region."""
    
    protein_id: str           # ID da proteína (qseqid)
    hsps: List[BlastHit]      # Lista de HSPs
    total_score: int          # Soma dos scores
    total_length: int         # Soma dos comprimentos das sequências
    score_density: float      # Densidade de score
    
    @classmethod
    def from_hsps(cls, protein_id: str, hsps: List[BlastHit]) -> 'ProteinHitGroup':
        """
        Create a ProteinHitGroup from a list of HSPs.
        
        Follows expert specification:
        Score Density = Σ(Scores of HSPs) / Σ(Alignment Lengths of HSPs)
        
        Where:
        - Numerator: Sum of raw scores from each HSP
        - Denominator: Sum of actual alignment lengths (hsp.length)
          NOT query coverage (qend-qstart) which would include gaps between HSPs
        
        Args:
            protein_id: Protein ID
            hsps: List of HSPs for this protein
        
        Returns:
            ProteinHitGroup object with calculated density
        """
        total_score = sum(hsp.score for hsp in hsps)
        # Use actual alignment length (number of aligned residues)
        # NOT query coverage which includes gaps between HSPs
        total_length = sum(hsp.length for hsp in hsps)
        
        # Calculate density: score / alignment_length
        score_density = total_score / total_length if total_length > 0 else 0.0
        
        return cls(
            protein_id=protein_id,
            hsps=hsps,
            total_score=total_score,
            total_length=total_length,
            score_density=score_density
        )
    
    def to_dict(self) -> Dict:
        return {
            'protein_id': self.protein_id,
            'hsps': [hsp.to_dict() for hsp in self.hsps],
            'total_score': self.total_score,
            'total_length': self.total_length,
            'score_density': self.score_density,
        }

    @classmethod
    def from_dict(cls, d: Dict) -> 'ProteinHitGroup':
        return cls(
            protein_id=d['protein_id'],
            hsps=[BlastHit.from_dict(h) for h in d['hsps']],
            total_score=d['total_score'],
            total_length=d['total_length'],
            score_density=d['score_density'],
        )


@dataclass
class RegionAnnotation:
    """Complete annotation of a genomic region."""
    
    region: GenomicRegion           # Região genômica fundida
    best_protein: ProteinHitGroup   # Proteína com maior densidade de score
    all_proteins: List[ProteinHitGroup]  # Todas as proteínas que mapearam

    def to_dict(self) -> Dict:
        """Serialização completa para JSON."""
        return {
            'region': self.region.to_dict(),
            'best_protein': self.best_protein.to_dict(),
            'all_proteins': [p.to_dict() for p in self.all_proteins],
        }
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'RegionAnnotation':
        """Deserializa um RegionAnnotation a partir de um dicionário."""
        return cls(
            region=GenomicRegion.from_dict(d['region']),
            best_protein=ProteinHitGroup.from_dict(d['best_protein']),
            all_proteins=[ProteinHitGroup.from_dict(p) for p in d['all_proteins']],
        )
    
    # def to_dict(self) -> Dict:
    #     """Convert annotation to dictionary."""
    #     return {
    #         'chromosome': self.region.chromosome,
    #         'start': self.region.start,
    #         'end': self.region.end,
    #         'length': self.region.length,
    #         'strand': self.region.strand,
    #         'best_protein_id': self.best_protein.protein_id,
    #         'best_protein_score_density': self.best_protein.score_density,
    #         'best_protein_total_score': self.best_protein.total_score,
    #         'best_protein_total_length': self.best_protein.total_length,
    #         'best_protein_num_hsps': len(self.best_protein.hsps),
    #         'num_competing_proteins': len(self.all_proteins)
    #     }


def calculate_score_density(hsps: List[BlastHit]) -> float:
    """
    Calculate score density for a set of HSPs.
    
    Follows expert specification:
    Density = Σ(Score_HSPs) / Σ(Alignment_Length_HSPs)
    
    Where:
    - Score: Raw BLAST/SSEARCH score for each HSP
    - Alignment Length: Actual number of aligned residues (hsp.length)
      NOT query coverage which would include unaligned gaps
    
    This measures alignment quality per aligned residue, answering:
    "In the regions where alignment occurred, how strong was that alignment?"
    
    Example:
      HSP1: score=150, length=50 aa
      HSP2: score=120, length=40 aa
      Density = (150+120)/(50+40) = 270/90 = 3.0
    
    Args:
        hsps: List of HSPs (High Scoring Pairs)
    
    Returns:
        Score density (score per aligned amino acid)
    """
    if not hsps:
        return 0.0
    
    total_score = sum(hsp.score for hsp in hsps)
    # Use actual alignment length (number of aligned residues)
    total_length = sum(hsp.length for hsp in hsps)
    
    density = total_score / total_length if total_length > 0 else 0.0
    
    logger.debug(
        f"Density calculated: {density:.2f} "
        f"(score={total_score}, alignment_length={total_length}, num_hsps={len(hsps)})"
    )
    
    return density


def group_hsps_by_protein(hsps: List[BlastHit]) -> Dict[str, List[BlastHit]]:
    """
    Group HSPs by query protein.
    
    Args:
        hsps: List of HSPs
    
    Returns:
        Dictionary mapping protein_id to list of its HSPs
    """
    protein_groups = defaultdict(list)
    
    for hsp in hsps:
        protein_groups[hsp.qseqid].append(hsp)
    
    return dict(protein_groups)


def select_best_protein_for_region(
    region: GenomicRegion,
    region_hsps: List[BlastHit]
) -> RegionAnnotation:
    """
    Select protein with highest score density for a region.
    
    Args:
        region: Merged genomic region
        region_hsps: List of HSPs that map to this region
    
    Returns:
        Region annotation with best protein selected
    """
    logger.debug(
        f"Selecting best protein for region "
        f"{region.chromosome}:{region.start}-{region.end}"
    )
    
    # Group HSPs by protein
    protein_groups = group_hsps_by_protein(region_hsps)
    
    # Calculate density for each protein
    protein_hit_groups = []
    for protein_id, hsps in protein_groups.items():
        hit_group = ProteinHitGroup.from_hsps(protein_id, hsps)
        protein_hit_groups.append(hit_group)
    
    # Select protein with highest density ("King of the Hill" algorithm)
    # Primary criterion: Highest Score Density
    # Tie-breaker: Highest Total Score (if densities are equal)
    best_protein = max(
        protein_hit_groups,
        key=lambda x: (x.score_density, x.total_score)
    )
    
    logger.debug(
        f"Best protein for region {region.chromosome}:{region.start}-{region.end}: "
        f"{best_protein.protein_id} (density={best_protein.score_density:.2f}, "
        f"total_score={best_protein.total_score})"
    )
    
    # Update query_ids list in region
    region.query_ids = [best_protein.protein_id]
    
    return RegionAnnotation(
        region=region,
        best_protein=best_protein,
        all_proteins=protein_hit_groups
    )


def assign_hsps_to_regions(
    regions: List[GenomicRegion],
    all_hsps: List[BlastHit]
) -> Dict[Tuple[str, int, int], List[BlastHit]]:
    """
    Assign HSPs to merged genomic regions.
    
    Args:
        regions: List of merged genomic regions
        all_hsps: List of all HSPs
    
    Returns:
        Dictionary mapping (chromosome, start, end) to list of HSPs
    """
    logger.debug(f"Assigning {len(all_hsps)} HSPs to {len(regions)} regions")
    
    region_hsps_map = defaultdict(list)
    
    for hsp in all_hsps:
        # HSP coordinates
        hsp_chrom = hsp.sseqid
        hsp_start = min(hsp.sstart, hsp.send)
        hsp_end = max(hsp.sstart, hsp.send)
        
        # Find regions that overlap with this HSP
        for region in regions:
            if region.chromosome != hsp_chrom:
                continue
            
            # Check overlap
            if hsp_start <= region.end and hsp_end >= region.start:
                region_key = (region.chromosome, region.start, region.end)
                region_hsps_map[region_key].append(hsp)
    
    logger.debug(f"HSPs assigned to {len(region_hsps_map)} regions")
    return dict(region_hsps_map)

import json

def annotate_regions_with_best_proteins(
    regions: List[GenomicRegion],
    all_hsps: List[BlastHit],
    output_path: Path,
) -> List[RegionAnnotation]:
    """
    Annotate each genomic region with the best protein (highest score density) 
    and saves them in a JSON file.
    
    Args:
        regions: List of merged genomic regions
        all_hsps: List of all filtered HSPs
        output_path: Path to output file in JSON format    
    Returns:
        List of region annotations with best proteins
    """
    logger.debug(f"Annotating {len(regions)} regions with best proteins...")
    
    # Assign HSPs to regions
    region_hsps_map = assign_hsps_to_regions(regions, all_hsps)
    
    # For each region, select best protein
    annotations = []
    for region in regions:
        region_key = (region.chromosome, region.start, region.end)
        region_hsps = region_hsps_map.get(region_key, [])
        
        if not region_hsps:
            logger.warning(
                f"No HSP found for region "
                f"{region.chromosome}:{region.start}-{region.end}"
            )
            continue
        
        annotation = select_best_protein_for_region(region, region_hsps)
        annotations.append(annotation)
    
    logger.debug(f"{len(annotations)} regions annotated")
    
    logger.debug(f"Saving {len(annotations)} annotations to: {output_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w') as f:
        for ann in annotations:

            logger.info(f"{ann}\n\n")
            f.write(json.dumps(ann.to_dict()) + '\n')  # JSON Lines: 1 objeto por linha

    logger.debug("Annotations saved successfully in JSONL")

    return annotations


import sys
from src.ssearch_realign import load_blast_hits_from_tsv
from src.bedtools_merge import parse_merged_bed
if __name__ == '__main__':
    genomic_region = parse_merged_bed(Path(sys.argv[1]), genome_id=sys.argv[4])
    blast_hits = load_blast_hits_from_tsv(Path(sys.argv[2]))
    annotate_regions_with_best_proteins(genomic_region, blast_hits, Path(sys.argv[3]))