import json
from pathlib import Path
from typing import List
from src.pseudogene_detection import PseudogeneAnnotation

def apply_coverage(pseudogene_annotations: List[PseudogeneAnnotation], min_coverage: float, output_path: Path) -> List[PseudogeneAnnotation]:
    filtered_annotations = []
        
    for ann in pseudogene_annotations:
        hsps = ann.region_annotation.best_protein.hsps
        if hsps:
            min_qstart = min(hsp.qstart for hsp in hsps)
            max_qend = max(hsp.qend for hsp in hsps)
            coverage_len = max_qend - min_qstart + 1
            ref_len = hsps[0].qlen
            coverage_ratio = coverage_len / ref_len if ref_len > 0 else 0
            
            if coverage_ratio >= min_coverage:
                filtered_annotations.append(ann)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w') as f:
        for ann in filtered_annotations:
            f.write(json.dumps(ann.to_json_dict()) + '\n')  # JSON Lines: 1 objeto por linha
    
    return filtered_annotations

def load_pseudogene_annotations_from_json(json_path: Path) -> List[PseudogeneAnnotation]:
    """Carrega PseudogeneAnnotations de um arquivo JSON Lines."""
    annotations = []
    with json_path.open('r') as f:
        for line in f:
            line = line.strip()
            if line:
                d = json.loads(line)
                annotations.append(PseudogeneAnnotation.from_dict(d))
    return annotations

import sys
if __name__ == "__main__":
    pseudogene_annotation = load_pseudogene_annotations_from_json(Path(sys.argv[1]))
    apply_coverage(pseudogene_annotation, float(sys.argv[2]), Path(sys.argv[3]))