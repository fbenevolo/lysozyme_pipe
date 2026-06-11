from pathlib import Path
from typing import List
from src.blast_filter import BlastHit
from src.ssearch_realign import SSearchAlignment

from src.ssearch_realign import load_blast_hits_from_tsv # convert TSV to List[BlastHit]
from src.filter_ssearch import load_alignments_from_tsv # convert TSV to Dict[str, SSearchAlignment]

def filter_hits_after_ssearch(filtered_hits: List[BlastHit], filtered_ssearch: dict[str, SSearchAlignment],
                              ssearch_hit_keys: List, filtered_hits_after_ssearch_output: Path) -> List[BlastHit]:
    """
    TODO generate docstring
        
    :param filtered_hits: Description
    :type filtered_hits: List[BlastHit]
    :param filtered_ssearch: Description
    :type filtered_ssearch: dict[str, SSearchAlignment]
    :param ssearch_hit_keys: Description
    :type ssearch_hit_keys: List
    :param filtered_hits_after_ssearch_output: Description
    :type filtered_hits_after_ssearch_output: Path
    :return: Description
    :rtype: List[BlastHit]
    """
    filtered_hits_after_ssearch = []
    
    for hit in filtered_hits:
        hit_key = f"{hit.qseqid}_{hit.sseqid}_{hit.sstart}_{hit.send}"
        if hit_key in ssearch_hit_keys:
            # Update BlastHit with SSEARCH scores
            ssearch_aln = filtered_ssearch[hit_key]
            
            # Update scores (critical for score density calculation)
            hit.score = int(ssearch_aln.bit_score)  # Use SSEARCH bit score as raw score
            hit.bitscore = ssearch_aln.bit_score
            hit.evalue = ssearch_aln.evalue
            
            # Update alignment details
            hit.length = ssearch_aln.alignment_length
            hit.pident = ssearch_aln.identity
            hit.mismatch = ssearch_aln.mismatches
            hit.gapopen = ssearch_aln.gap_opens
            
            filtered_hits_after_ssearch.append(hit)

    with open(filtered_hits_after_ssearch_output, 'w') as f:
        header = [
            "genome_id", "qseqid", "qlen", "sseqid", "slen", "qstart", "qend", 
            "sstart", "send", "qseq", "sseq", "evalue", "bitscore", "score", 
            "length", "pident", "nident", "mismatch", "positive", "gapopen", 
            "gaps", "ppos", "sframe", "sstrand", "qcovs", "qcovhsp"
        ]
        f.write('\t'.join(header) + '\n')
        
        for hit in filtered_hits_after_ssearch:
            line = [
                str(hit.genome_id),
                str(hit.qseqid),
                str(hit.qlen),
                str(hit.sseqid),
                str(hit.slen),
                str(hit.qstart),
                str(hit.qend),
                str(hit.sstart),
                str(hit.send),
                str(hit.qseq),
                str(hit.sseq),
                f"{hit.evalue:.2e}",
                f"{hit.bitscore:.1f}",
                str(hit.score),
                str(hit.length),
                f"{hit.pident:.2f}",
                str(hit.nident),
                str(hit.mismatch),
                str(hit.positive),
                str(hit.gapopen),
                str(hit.gaps),
                f"{hit.ppos:.2f}" if isinstance(hit.ppos, (int, float)) else str(hit.ppos),
                str(hit.sframe),
                str(hit.sstrand),
                str(hit.qcovs),
                str(hit.qcovhsp)
            ]
            f.write('\t'.join(line) + '\n')


    return filtered_hits_after_ssearch


import sys
if __name__ == '__main__':
    blast_hit_list = load_blast_hits_from_tsv(sys.argv[1])
    ssearch_aligment_to_hit = load_alignments_from_tsv(sys.argv[2])
    ssearch_hit_keys = ssearch_aligment_to_hit.keys()
    filter_hits_after_ssearch(blast_hit_list, ssearch_aligment_to_hit, ssearch_hit_keys, sys.argv[3])