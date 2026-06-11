"""
Module for exporting annotations in GFF3 format.
Generates GFF3 files compatible with genomic visualization tools (IGV, JBrowse, Artemis).
"""

import logging
from pathlib import Path
from typing import List
from datetime import datetime

from src.pseudogene_detection import PseudogeneAnnotation


logger = logging.getLogger(__name__)


def export_to_gff3(
    pseudogene_annotations: List[PseudogeneAnnotation],
    genome_id: str,
    output_path: Path,
    source: str = "lysozyme_pipeline"
) -> None:
    """
    Export pseudogene annotations to GFF3 format.
    
    GFF3 format:
    seqid source type start end score strand phase attributes
    
    Args:
        pseudogene_annotations: List of pseudogene annotations
        genome_id: Genome identifier
        output_path: Path to output GFF3 file
        source: Tool/pipeline name (column 2)
    """
    logger.debug(f"Exporting {len(pseudogene_annotations)} annotations to GFF3: {output_path}")
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    
    with open(output_path, 'w') as f:
        # Write GFF3 header
        f.write("##gff-version 3\n")
        f.write(f"##date {datetime.now().strftime('%Y-%m-%d')}\n")
        f.write(f"##source lysozyme_pipeline v1.0\n")
        f.write(f"##genome-id {genome_id}\n\n")
        
        # Write each annotation
        for i, annotation in enumerate(pseudogene_annotations, 1):
            # Determine feature type
            feature_type = "pseudogene" if annotation.is_pseudogene else "CDS"
            
            # Determine strand (assume + if no information)
            strand = "+"
            
            # Determine score (use best protein score_density)
            score_density = annotation.region_annotation.best_protein.score_density
            score = f"{score_density:.2f}" if score_density > 0 else "."
            
            # Build attributes column (9)
            attributes = _build_gff3_attributes(annotation, genome_id, i)
            
            # GFF3 line
            gff3_line = "\t".join([
                annotation.region_annotation.region.chromosome,  # seqid
                source,                                   # source
                feature_type,                             # type
                str(annotation.region_annotation.region.start),  # start
                str(annotation.region_annotation.region.end),    # end
                score,                                    # score
                strand,                                   # strand
                ".",                                      # phase (not applicable)
                attributes                                # attributes
            ])
            
            f.write(gff3_line + "\n")
    
    logger.debug(f"GFF3 file saved: {output_path}")


def _build_gff3_attributes(
    annotation: PseudogeneAnnotation,
    genome_id: str,
    index: int
) -> str:
    """
    Build attributes column (column 9) for GFF3.
    
    Format: key1=value1;key2=value2;...
    
    Args:
        annotation: Pseudogene annotation
        genome_id: Genome identifier
        index: Annotation index (for unique ID)
    
    Returns:
        String with formatted attributes
    """
    attributes = []
    
    # Unique feature ID
    feature_id = f"lysozyme_{genome_id}_{index}"
    attributes.append(f"ID={feature_id}")
    
    # Parent (genome)
    attributes.append(f"Parent={genome_id}")
    
    # Functional status
    status = "Pseudogene" if annotation.is_pseudogene else "Functional"
    attributes.append(f"Status={status}")
    
    # Small ORF status
    if annotation.is_small_orf:
        attributes.append("Small_ORF=True")
    
    # Reference protein
    best_protein_id = annotation.region_annotation.best_protein.protein_id
    if best_protein_id:
        attributes.append(f"Ref_Protein={best_protein_id}")
    
    # Detected mutation types
    mutations = []
    if annotation.disablements.non_synonymous_substitutions > 0:
        mutations.append(f"NonSynSubst:{annotation.disablements.non_synonymous_substitutions}")
    if annotation.disablements.in_frame_indels > 0:
        mutations.append(f"InFrameIndels:{annotation.disablements.in_frame_indels}")
    if annotation.disablements.frameshifts > 0:
        mutations.append(f"Frameshifts:{annotation.disablements.frameshifts}")
    if annotation.disablements.premature_stop_codons > 0:
        mutations.append(f"PrematureStops:{annotation.disablements.premature_stop_codons}")
    if annotation.disablements.missing_start_codon > 0:
        mutations.append(f"MissingStart:{annotation.disablements.missing_start_codon}")
    if annotation.disablements.missing_stop_codon > 0:
        mutations.append(f"MissingStop:{annotation.disablements.missing_stop_codon}")
    
    if mutations:
        attributes.append(f"Mutations={','.join(mutations)}")
    
    # Total inactivating mutations
    attributes.append(f"Total_Disablements={annotation.disablements.total_disablements}")
    
    # Score density
    score_density = annotation.region_annotation.best_protein.score_density
    attributes.append(f"Score_Density={score_density:.2f}")
    
    # Descriptive name
    name = f"Lysozyme_homolog_{index}"
    if annotation.is_pseudogene:
        name += "_pseudogene"
    attributes.append(f"Name={name}")
    
    return ";".join(attributes)


def export_batch_to_gff3(
    annotations_by_genome: dict,
    output_dir: Path,
    source: str = "lysozyme_pipeline"
) -> None:
    """
    Export annotations from multiple genomes to individual GFF3 files.
    
    Args:
        annotations_by_genome: Dict mapping genome_id to list of annotations
        output_dir: Base output directory
        source: Tool/pipeline name
    """
    logger.debug(f"Exporting annotations from {len(annotations_by_genome)} genomes to GFF3")
    
    for genome_id, annotations in annotations_by_genome.items():
        genome_output_path = output_dir / f"{genome_id}_annotations.gff3"
        export_to_gff3(annotations, genome_id, genome_output_path, source)
    
    logger.debug(f"Batch export complete. Files in: {output_dir}")



import sys
from src.save_pseudogene_annotations import load_pseudogene_annotations_from_json
if __name__ == "__main__":
    pseudogene_annotations = load_pseudogene_annotations_from_json(Path(sys.argv[1]))
    export_to_gff3(pseudogene_annotations, sys.argv[2], Path(sys.argv[3]))