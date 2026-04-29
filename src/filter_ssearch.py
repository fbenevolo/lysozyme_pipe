import logging
from ssearch_realign import SSearchAlignment
from pathlib import Path
from typing import Dict

logger = logging.getLogger(__name__)


def filter_ssearch_by_evalue(ssearch_alignments: dict[str, SSearchAlignment],
                             ssearch_output: Path) -> dict[str, SSearchAlignment]:
    evalue_threshold = 1e-7
    logger.debug(f"Filtering SSEARCH by E-value <= {evalue_threshold:.0e}")

    
    filtered_ssearch = {key: aln for key, aln in ssearch_alignments.items()
                       if aln.evalue <= evalue_threshold}
    
    with open(ssearch_output, 'w') as f:
        f.write("hit_key\tquery_id\tsubject_id\tidentity\tevalue\tbit_score\n")
        for key, aln in filtered_ssearch.items():
            f.write(f"{key}\t{aln.query_id}\t{aln.subject_id}\t"
                   f"{aln.identity:.2f}\t{aln.evalue:.2e}\t{aln.bit_score:.2f}\n")
    
    return filtered_ssearch


import pandas as pd
import sys

def load_alignments_from_tsv(tsv_path: Path) -> Dict[str, SSearchAlignment]:
    df = pd.read_csv(tsv_path, sep='\t')
    alignments = {
        str(row['hit_key']): SSearchAlignment(
            **{k: v for k, v in row.items() if k != 'hit_key'}
        )
        for row in df.to_dict('records')
    }
    return alignments

if __name__ == '__main__':
    tsv_file = sys.argv[1]
    alignments = load_alignments_from_tsv(tsv_file)
    filter_ssearch_by_evalue(alignments, Path(sys.argv[2]))