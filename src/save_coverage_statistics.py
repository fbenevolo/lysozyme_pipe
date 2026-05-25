import logging
from pathlib import Path
from typing import List
from pseudogene_detection import PseudogeneAnnotation, calculate_subject_coverage_nt

logger = logging.getLogger(__name__)

def save_coverage_statistics(
    annotations: List[PseudogeneAnnotation],
    output_path
) -> None:
    """
    Save detailed coverage statistics for further analysis.
    
    Creates TSV with coverage ratio, alignment details, and classification
    for discussion about size-based validation approaches.
    
    Args:
        annotations: List of pseudogene annotations
        output_path: Path to coverage statistics output file
    """
    import pandas as pd
    
    logger.debug(f"Saving coverage statistics to: {output_path}")

    coverage_data = []
    for ann in annotations:
        hsps = ann.region_annotation.best_protein.hsps
        if hsps:
            min_qstart = min(hsp.qstart for hsp in hsps)
            max_qend = max(hsp.qend for hsp in hsps)
            coverage_len = max_qend - min_qstart + 1
            ref_len = hsps[0].qlen
            coverage_ratio = coverage_len / ref_len if ref_len > 0 else 0
            
            # Calculate genomic length
            genomic_len_nt = ann.region_annotation.region.length
            genomic_len_aa = genomic_len_nt / 3
            
            # Calculate subject coverage (nt)
            alignment_coverage_nt = calculate_subject_coverage_nt(hsps)
            
            # Calculate reference coverage in nt (alignment span on reference * 3)
            reference_coverage_nt = coverage_len * 3
            
            # Calculate alignment/genomic ratio (User requested metric)
            # (reference_coverage_nt) / genomic_region_nt
            alignment_genomic_ratio = reference_coverage_nt / genomic_len_nt if genomic_len_nt > 0 else 0
            
            # Determine classification
            if ann.is_pseudogene:
                if ann.is_small_orf:
                    classification = 'Pseudogene (Small ORF)'
                else:
                    classification = 'Pseudogene (Detected)'
            else:
                if ann.is_small_orf:
                    classification = 'Functional (Small ORF)'
                else:
                    classification = 'Functional (Possible Gene)'
            
            coverage_data.append({
                'chromosome': ann.region_annotation.region.chromosome,
                'start': ann.region_annotation.region.start,
                'end': ann.region_annotation.region.end,
                'strand': ann.region_annotation.region.strand,
                'protein_id': ann.region_annotation.best_protein.protein_id,
                'is_pseudogene': ann.is_pseudogene,
                'is_small_orf': ann.is_small_orf,
                'classification': classification,
                'alignment_coverage_aa': coverage_len,
                'alignment_coverage_nt': alignment_coverage_nt,
                'reference_coverage_nt': reference_coverage_nt,
                'alignment_genomic_ratio': round(alignment_genomic_ratio, 3),
                'reference_length_aa': ref_len,
                'coverage_ratio': round(coverage_ratio, 3),
                'genomic_region_nt': genomic_len_nt,
                'genomic_region_aa': round(genomic_len_aa, 1),
                'num_hsps': len(hsps),
                'qstart_min': min_qstart,
                'qend_max': max_qend,
                'missing_start_codon': ann.disablements.missing_start_codon,
                'missing_stop_codon': ann.disablements.missing_stop_codon,
                'premature_stops': ann.disablements.premature_stop_codons,
                'frameshifts': ann.disablements.frameshifts
            })
    
    df = pd.DataFrame(coverage_data)
    
    # Sort by coverage ratio (ascending) to highlight problematic cases
    df = df.sort_values('coverage_ratio', ascending=True)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)
    
    logger.debug(f"Coverage statistics saved: {len(coverage_data)} regions")


import sys
from apply_coverage_filter import load_pseudogene_annotations_from_json
if __name__ == "__main__":
    pseudogene_annotations = load_pseudogene_annotations_from_json(Path(sys.argv[1]))
    save_coverage_statistics(pseudogene_annotations, Path(sys.argv[2]))