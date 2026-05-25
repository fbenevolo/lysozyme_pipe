from pathlib import Path
from typing import List

from pseudogene_detection import PseudogeneAnnotation


def generate_summary_report(annotations: List[PseudogeneAnnotation], output_path: Path, min_coverage: float = 0.8) -> str:
    """
    Gera relatório resumido das anotações.
    
    Args:
        annotations: Lista de anotações
        min_coverage: Limite de cobertura para relatório (default: 0.8)
    
    Returns:
        String com o relatório formatado
    """
    total_regions = len(annotations)
    
    # Classification Breakdown
    functional_anns = [ann for ann in annotations if not ann.is_pseudogene]
    pseudogene_anns = [ann for ann in annotations if ann.is_pseudogene]
    
    num_functional = len(functional_anns)
    num_pseudogenes = len(pseudogene_anns)
    
    # Functional Sub-categories
    func_possible_genes = sum(1 for ann in functional_anns if not ann.is_small_orf)
    func_small_orfs = sum(1 for ann in functional_anns if ann.is_small_orf)
    
    # Pseudogene Sub-categories
    pseudo_detected = sum(1 for ann in pseudogene_anns if not ann.is_small_orf)
    pseudo_small_orfs = sum(1 for ann in pseudogene_anns if ann.is_small_orf)
    
    # Estatísticas de mutações
    total_substitutions = sum(ann.disablements.non_synonymous_substitutions for ann in annotations)
    total_indels = sum(ann.disablements.in_frame_indels for ann in annotations)
    total_frameshifts = sum(ann.disablements.frameshifts for ann in annotations)
    total_missing_start = sum(ann.disablements.missing_start_codon for ann in annotations)
    total_missing_stop = sum(ann.disablements.missing_stop_codon for ann in annotations)
    total_premature_stops = sum(ann.disablements.premature_stop_codons for ann in annotations)
    
    # NOVA SEÇÃO: Análise de cobertura e tamanho
    coverage_stats = []
    for ann in annotations:
            
        hsps = ann.region_annotation.best_protein.hsps
        if hsps:
            min_qstart = min(hsp.qstart for hsp in hsps)
            max_qend = max(hsp.qend for hsp in hsps)
            coverage_len = max_qend - min_qstart + 1
            ref_len = hsps[0].qlen
            coverage_ratio = coverage_len / ref_len if ref_len > 0 else 0
            
            coverage_stats.append({
                'coverage_len': coverage_len,
                'ref_len': ref_len,
                'ratio': coverage_ratio,
                'is_pseudogene': ann.is_pseudogene,
                'is_small_orf': ann.is_small_orf
            })
    
    # Helper function for stats
    def calc_stats(ratios):
        if not ratios:
            return 0, 0, 0, 0
        avg = sum(ratios) / len(ratios)
        mn = min(ratios)
        mx = max(ratios)
        below_threshold = sum(1 for r in ratios if r < min_coverage)
        return avg, mn, mx, below_threshold
    
    # DEBUG
    # print(f"DEBUG: min_coverage={min_coverage}")

    # 1. Functional - Possible Genes
    func_possible_ratios = [s['ratio'] for s in coverage_stats if not s['is_pseudogene'] and not s['is_small_orf']]
    fp_avg, fp_min, fp_max, fp_below = calc_stats(func_possible_ratios)
    fp_count = len(func_possible_ratios)

    # 2. Functional - Small ORFs
    func_small_ratios = [s['ratio'] for s in coverage_stats if not s['is_pseudogene'] and s['is_small_orf']]
    fs_avg, fs_min, fs_max, fs_below = calc_stats(func_small_ratios)
    fs_count = len(func_small_ratios)

    # 3. Pseudogenes - Detected
    pseudo_detected_ratios = [s['ratio'] for s in coverage_stats if s['is_pseudogene'] and not s['is_small_orf']]
    pd_avg, pd_min, pd_max, pd_below = calc_stats(pseudo_detected_ratios)
    pd_count = len(pseudo_detected_ratios)

    # 4. Pseudogenes - Small ORFs
    pseudo_small_ratios = [s['ratio'] for s in coverage_stats if s['is_pseudogene'] and s['is_small_orf']]
    ps_avg, ps_min, ps_max, ps_below = calc_stats(pseudo_small_ratios)
    ps_count = len(pseudo_small_ratios)
    
    report = f"""
╭──────────────────────────────────────────────────────────────────╮
│          LYSOZYME PSEUDOGENE ANNOTATION REPORT                   │
╰──────────────────────────────────────────────────────────────────╯

GENERAL SUMMARY:
  Total regions analyzed:          {total_regions}
  
  Classification:
    Functional Genes:                {num_functional} ({100*num_functional/total_regions if total_regions else 0:.1f}%)
      - Possible Genes:              {func_possible_genes} ({100*func_possible_genes/num_functional if num_functional else 0:.1f}%)
      - Small ORFs:                  {func_small_orfs} ({100*func_small_orfs/num_functional if num_functional else 0:.1f}%)
      
    Pseudogenes:                     {num_pseudogenes} ({100*num_pseudogenes/total_regions if total_regions else 0:.1f}%)
      - Detected Pseudogenes:        {pseudo_detected} ({100*pseudo_detected/num_pseudogenes if num_pseudogenes else 0:.1f}%)
      - Small ORFs:                  {pseudo_small_orfs} ({100*pseudo_small_orfs/num_pseudogenes if num_pseudogenes else 0:.1f}%)

MUTATION STATISTICS:
  Non-synonymous substitutions:    {total_substitutions}
  In-frame indels:                 {total_indels}
  Frameshifts:                     {total_frameshifts}
  Missing start codon:             {total_missing_start}
  Missing stop codon:              {total_missing_stop}
  Premature stop codons:           {total_premature_stops}
  
  Total inactivating mutations:    {total_substitutions + total_indels + total_frameshifts + total_missing_start + total_missing_stop + total_premature_stops}

REFERENCE PROTEIN COVERAGE ANALYSIS:
  (Ratio = Alignment Coverage / Reference Size)
  
  1. Functional - Possible Genes ({fp_count} regions):
    Mean coverage ratio:             {fp_avg:.2f}
    Minimum ratio:                   {fp_min:.2f}
    Maximum ratio:                   {fp_max:.2f}
    Regions with coverage <{int(min_coverage*100)}%:      {fp_below} ({100*fp_below/fp_count if fp_count else 0:.1f}%)

  2. Functional - Small ORFs ({fs_count} regions):
    Mean coverage ratio:             {fs_avg:.2f}
    Minimum ratio:                   {fs_min:.2f}
    Maximum ratio:                   {fs_max:.2f}
    Regions with coverage <{int(min_coverage*100)}%:      {fs_below} ({100*fs_below/fs_count if fs_count else 0:.1f}%)
  
  3. Pseudogenes - Detected ({pd_count} regions):
    Mean coverage ratio:             {pd_avg:.2f}
    Minimum ratio:                   {pd_min:.2f}
    Maximum ratio:                   {pd_max:.2f}
    Regions with coverage <{int(min_coverage*100)}%:      {pd_below} ({100*pd_below/pd_count if pd_count else 0:.1f}%)

  4. Pseudogenes - Small ORFs ({ps_count} regions):
    Mean coverage ratio:             {ps_avg:.2f}
    Minimum ratio:                   {ps_min:.2f}
    Maximum ratio:                   {ps_max:.2f}
    Regions with coverage <{int(min_coverage*100)}%:      {ps_below} ({100*ps_below/ps_count if ps_count else 0:.1f}%)


"""
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, 'w') as f:
        f.write(report)

    
    return report


import sys
from apply_coverage_filter import load_pseudogene_annotations_from_json
if __name__ == "__main__":
    pseudogene_annotations = load_pseudogene_annotations_from_json(Path(sys.argv[1]))
    generate_summary_report(pseudogene_annotations, Path(sys.argv[2]), float(sys.argv[3]))