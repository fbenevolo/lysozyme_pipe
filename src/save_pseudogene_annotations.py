import json
import pandas as pd
from pathlib import Path
from typing import List
from src.pseudogene_detection import PseudogeneAnnotation


def save_pseudogene_annotations(
    annotations: List[PseudogeneAnnotation],
    output_path
) -> None:
    """
    Salva anotações de pseudogenes em arquivo TSV.
    
    Args:
        annotations: Lista de anotações de pseudogenes
        output_path: Caminho para o arquivo de saída
    """
    
    data = [ann.to_tsv_dict() for ann in annotations]
    df = pd.DataFrame(data)
    
    output_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(output_path, sep='\t', index=False)

def load_pseudogene_annotations_from_json(annotations_path: Path) -> List[PseudogeneAnnotation]:
    annotations = []
    with annotations_path.open('r') as f:
        for line in f:
            line = line.strip()
            if line:
                d = json.loads(line)
                annotations.append(PseudogeneAnnotation.from_dict(d))
    return annotations


import sys
if __name__ == "__main__":
    pseudogene_annotations = load_pseudogene_annotations_from_json(Path(sys.argv[1]))
    save_pseudogene_annotations(pseudogene_annotations, Path(sys.argv[2]))