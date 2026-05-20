import json
from pathlib import Path
from typing import List
from pseudogene_detection import PseudogeneAnnotation


def apply_final_identity(pseudogene_annotations: List[PseudogeneAnnotation], final_min_identity: float,
                         output_path: Path) -> List[PseudogeneAnnotation]:
    threshold_pct = final_min_identity * 100 if final_min_identity <= 1.0 else final_min_identity
    
    filtered_annotations = []
    for ann in pseudogene_annotations:
        if max(hsp.pident for hsp in ann.region_annotation.best_protein.hsps) >= threshold_pct:
            filtered_annotations.append(ann)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
         for ann in filtered_annotations:
            f.write(json.dumps(ann.to_dict()) + '\n')  # JSON Lines: 1 objeto por linha

    return filtered_annotations

import sys
from apply_coverage_filter import load_pseudogene_annotations_from_json
if __name__ == "__main__":
    pseudogene_annotation = load_pseudogene_annotations_from_json(Path(sys.argv[1]))
    apply_final_identity(pseudogene_annotation, float(sys.argv[2]), Path(sys.argv[3]))