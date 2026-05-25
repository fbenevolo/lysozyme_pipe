"""
Módulo para detecção de mutações e caracterização de pseudogenes.
Implementa a etapa 6 do pipeline: Detecção de Mutação e Pseudogenes.
"""

import json
import logging
from bedtools_merge import GenomicRegion
import pandas as pd
from typing import List, Set, Dict, Optional
from dataclasses import dataclass
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq

from src.config import (
    DisablementCounts,
    GAP_CHAR,
    STOP_CODON_CHAR,
    START_CODON_AA,
    MIN_LENGTH_RATIO,
    MAX_LENGTH_RATIO,
    MIN_REGION_SIZE_NT,
    MIN_ALIGNMENT_SIZE_AA
)
from src.blast_filter import BlastHit
from src.score_density import RegionAnnotation, ProteinHitGroup


logger = logging.getLogger(__name__)


@dataclass
class PseudogeneAnnotation:
    """
    Anotação completa de um pseudogene potencial.
    
    Estrutura de acesso aos atributos:
    - annotation.is_pseudogene
    - annotation.region_annotation.region.chromosome
    - annotation.region_annotation.region.start
    - annotation.region_annotation.region.end
    - annotation.region_annotation.region.strand
    - annotation.region_annotation.best_protein.protein_id
    - annotation.region_annotation.best_protein.score_density
    - annotation.disablements.non_synonymous_substitutions
    - annotation.disablements.in_frame_indels
    - annotation.disablements.frameshifts
    - annotation.disablements.premature_stop_codons
    - annotation.disablements.missing_start_codon
    - annotation.disablements.missing_stop_codon
    - annotation.disablements.total_disablements (propriedade)
    """
    
    region_annotation: RegionAnnotation  # Anotação da região
    disablements: DisablementCounts      # Contadores de mutações
    is_pseudogene: bool                  # Se é classificado como pseudogene
    is_small_orf: bool = False           # If region is small (<300 nt / 100 aa)
    
    def to_json_dict(self):
        """Converte a anotação em dicionário."""
        result = self.region_annotation.to_json_dict()
        result.update(self.disablements.to_dict())
        result['is_pseudogene'] = self.is_pseudogene
        result['is_small_orf'] = self.is_small_orf
        
        # Add concatenated sequences for visualization
        hsps = self.region_annotation.best_protein.hsps
        if hsps:
            # Sort HSPs by query start to ensure correct order
            sorted_hsps = sorted(hsps, key=lambda h: h.qstart)
            result['qseq'] = "".join(h.qseq for h in sorted_hsps)
            result['sseq'] = "".join(h.sseq for h in sorted_hsps)
        else:
            result['qseq'] = ""
            result['sseq'] = ""
            
        return result
    
    def to_tsv_dict(self):
        result = self.region_annotation.to_tsv_dict()
        result.update(self.disablements.to_dict())
        result['is_pseudogene'] = self.is_pseudogene
        result['is_small_orf'] = self.is_small_orf
        
        # Add concatenated sequences for visualization
        hsps = self.region_annotation.best_protein.hsps
        if hsps:
            # Sort HSPs by query start to ensure correct order
            sorted_hsps = sorted(hsps, key=lambda h: h.qstart)
            result['qseq'] = "".join(h.qseq for h in sorted_hsps)
            result['sseq'] = "".join(h.sseq for h in sorted_hsps)
        else:
            result['qseq'] = ""
            result['sseq'] = ""
            
        return result
    
    # def to_dict(self) -> Dict:
    #     """Serialização completa para JSON — preserva estrutura aninhada."""
    #     return {
    #         'region_annotation': self.region_annotation.to_json_dict(),
    #         'disablements': self.disablements.to_dict(),
    #         'is_pseudogene': self.is_pseudogene,
    #         'is_small_orf': self.is_small_orf,
    #     }
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'PseudogeneAnnotation':
        # Reconstrói RegionAnnotation a partir da estrutura plana
        region = GenomicRegion(
            genome_id=d['region']['genome_id'],
            chromosome=d['region']['chromosome'],
            start=d['region']['start'],
            end=d['region']['end'],
            strand=d['region']['strand'],
            num_hsps=d['region']['num_hsps'],
            mean_score=d['region']['mean_score'],
            min_score=d['region']['min_score'],
            max_score=d['region']['max_score'],
            query_ids=d['region']['query_ids'].split(',') if isinstance(d['region']['query_ids'], str) else d['region']['query_ids'],
        )
        region_annotation = RegionAnnotation(
            region=region,
            best_protein=ProteinHitGroup.from_dict(d['best_protein']),
            all_proteins=[ProteinHitGroup.from_dict(p) for p in d['all_proteins']],
        )
        disablements = DisablementCounts(
            non_synonymous_substitutions=d['non_synonymous_substitutions'],
            in_frame_indels=d['in_frame_indels'],
            frameshifts=d['frameshifts'],
            missing_start_codon=d['missing_start_codon'],
            missing_stop_codon=d['missing_stop_codon'],
            premature_stop_codons=d['premature_stop_codons'],
            size_mismatch=d['size_mismatch'],
        )
        return cls(
            region_annotation=region_annotation,
            disablements=disablements,
            is_pseudogene=d['is_pseudogene'],
            is_small_orf=d.get('is_small_orf', False),
        )


def count_non_synonymous_substitutions(query_seq: str, subject_seq: str) -> int:
    """
    Conta substituições não-sinônimas entre sequências alinhadas.
    
    Conta aminoácidos não idênticos na mesma posição do alinhamento,
    excluindo gaps. IMPORTANTE: Nem toda substituição indica pseudogene,
    apenas contamos para análise.
    
    Args:
        query_seq: Sequência query alinhada (proteína funcional)
        subject_seq: Sequência subject alinhada (genômica traduzida)
    
    Returns:
        Número de substituições não-sinônimas
    """
    if len(query_seq) != len(subject_seq):
        logger.warning(
            f"Sequências com comprimentos diferentes: "
            f"query={len(query_seq)}, subject={len(subject_seq)}"
        )
    
    count = 0
    min_length = min(len(query_seq), len(subject_seq))
    
    for i in range(min_length):
        q_aa = query_seq[i]
        s_aa = subject_seq[i]
        
        # Ignora posições com gaps
        if q_aa == GAP_CHAR or s_aa == GAP_CHAR:
            continue
        
        # Ignora stop codons na contagem de substituições
        if q_aa == STOP_CODON_CHAR or s_aa == STOP_CODON_CHAR:
            continue
        
        # Conta se forem diferentes
        if q_aa != s_aa:
            count += 1
    
    # Calcula percentual de diferença
    total_positions = min_length - query_seq.count(GAP_CHAR) - subject_seq.count(GAP_CHAR)
    if total_positions > 0:
        percent_diff = (count / total_positions) * 100
        logger.debug(f"Substituições: {count}/{total_positions} ({percent_diff:.1f}%)")
    
    return count


def count_in_frame_indels(query_seq: str, subject_seq: str) -> int:
    """
    Conta indels in-frame (gaps individuais) presentes em ambas as sequências.
    
    Args:
        query_seq: Sequência query alinhada
        subject_seq: Sequência subject alinhada
    
    Returns:
        Número total de gaps individuais
    """
    query_gaps = query_seq.count(GAP_CHAR)
    subject_gaps = subject_seq.count(GAP_CHAR)
    
    total_gaps = query_gaps + subject_gaps
    
    return total_gaps


def count_frameshifts(hsps: List[BlastHit]) -> int:
    """
    Calcula o número de frameshifts baseado em frames distintos dos HSPs.
    
    Frameshifts = count(distinct_frames) - 1
    
    Args:
        hsps: Lista de HSPs da mesma sequência de referência
    
    Returns:
        Número de frameshifts
    """
    if not hsps:
        return 0
    
    # Coleta todos os frames únicos
    frames: Set[int] = set()
    for hsp in hsps:
        frames.add(hsp.sframe)
    
    # Número de frameshifts é o número de frames distintos menos 1
    num_frameshifts = len(frames) - 1
    
    logger.debug(
        f"Frames detectados: {sorted(frames)}, "
        f"Frameshifts: {num_frameshifts}"
    )
    
    return max(0, num_frameshifts)


def check_missing_start_codon(subject_seq: str) -> bool:
    """
    Verifica se o start codon (Metionina) está ausente.
    
    IMPORTANTE: tblastn traduz o DNA mas NÃO inclui o start codon (ATG->M) 
    automaticamente na sequência traduzida. Por isso, esta verificação é
    relevante apenas se o alinhamento começa no início real do gene.
    
    Args:
        subject_seq: Sequência subject alinhada (genômica traduzida)
    
    Returns:
        True se o start codon está ausente, False caso contrário
    """
    # Remove gaps do início
    trimmed_seq = subject_seq.lstrip(GAP_CHAR)
    
    if not trimmed_seq:
        return True
    
    # Para tblastn, a ausência de M no início pode ser normal se o alinhamento
    # não começa exatamente no start codon. Portanto, NÃO consideramos isso
    # como evidência de pseudogene por padrão.
    # Apenas reportamos se estiver completamente ausente E tiver outros indicadores
    first_aa = trimmed_seq[0]
    
    # Relaxa: só marca como problemático se não for M E a sequência for longa o suficiente
    # para esperar que comece no início real
    if len(trimmed_seq) < 50:  # Alinhamentos curtos podem não incluir o start
        return False
    
    is_missing = first_aa != START_CODON_AA
    
    if is_missing:
        logger.debug(f"Possível start codon ausente. Primeiro AA: {first_aa}")
    
    # NÃO marca como problema por padrão - retorna False
    return False


def check_missing_stop_codon(subject_seq: str) -> bool:
    """
    Verifica se o stop codon está ausente no final da sequência.
    
    IMPORTANTE: tblastn traduz até encontrar um stop codon, mas NÃO inclui
    o stop codon (*) na sequência traduzida resultante. Portanto, a ausência
    de * no final é NORMAL e NÃO indica pseudogene.
    
    Args:
        subject_seq: Sequência subject alinhada (genômica traduzida)
    
    Returns:
        True se o stop codon está ausente, False caso contrário
    """
    # Remove gaps do final
    trimmed_seq = subject_seq.rstrip(GAP_CHAR)
    
    if not trimmed_seq:
        return True
    
    # Para tblastn, o stop codon NÃO aparece na sequência traduzida normalmente.
    # A presença de * no final seria anormal. A ausência é esperada.
    # Portanto, SEMPRE retornamos False (não é problema)
    last_aa = trimmed_seq[-1]
    
    # Se houver * no final, isso seria estranho, mas não necessariamente problema
    if last_aa == STOP_CODON_CHAR:
        logger.debug(f"Stop codon presente no final (incomum mas não problemático): {last_aa}")
    
    # NÃO marca como problema - retorna sempre False
    return False


def count_premature_stop_codons(subject_seq: str) -> int:
    """
    Conta stop codons prematuros (internos) na sequência.
    
    Args:
        subject_seq: Sequência subject alinhada (genômica traduzida)
    
    Returns:
        Número de stop codons internos
    """
    # Remove gaps
    clean_seq = subject_seq.replace(GAP_CHAR, '')
    
    if not clean_seq:
        return 0
    
    # Remove o último caractere (que pode ser stop codon legítimo)
    internal_seq = clean_seq[:-1] if len(clean_seq) > 1 else ""
    
    # Conta stop codons internos
    count = internal_seq.count(STOP_CODON_CHAR)
    
    if count > 0:
        logger.debug(f"Stop codons prematuros detectados: {count}")
    
    return count


def check_length_viability(
    candidate_length: int,
    reference_length: int
) -> Dict[str, any]:
    """
    Apply the 20% Rule to validate gene functionality based on size.
    
    Validates whether candidate gene length is within acceptable range
    (80-120%) of reference protein length. Genes outside this range are
    likely truncated, elongated, or fusion artifacts.
    
    Args:
        candidate_length: Length of candidate sequence in amino acids
        reference_length: Length of reference protein in amino acids
    
    Returns:
        Dict with validation results:
            - candidate_len: Candidate length
            - ref_len: Reference length
            - ratio: Length ratio (candidate/reference)
            - status: "PASS" or "FAIL"
            - classification: Classification result
    """
    # Calculate thresholds (20% Rule)
    min_len = reference_length * MIN_LENGTH_RATIO
    max_len = reference_length * MAX_LENGTH_RATIO
    
    ratio = candidate_length / reference_length if reference_length > 0 else 0
    
    result = {
        "candidate_len": candidate_length,
        "ref_len": reference_length,
        "ratio": round(ratio, 2)
    }
    
    if min_len <= candidate_length <= max_len:
        result["status"] = "PASS"
        result["classification"] = "Size compatible"
    elif candidate_length < min_len:
        result["status"] = "FAIL"
        result["classification"] = "Truncated (size < 80%)"
    else:
        result["status"] = "FAIL"
        result["classification"] = "Elongated/Fusion (size > 120%)"
    
    return result


def check_start_stop_codons_genomic(
    genome_fasta_path: Path,
    chromosome: str,
    start: int,
    end: int,
    strand: str,
    reference_length_aa: int,
    check_start: bool = True,
    check_stop: bool = True
) -> Dict[str, bool]:
    """
    Check for start/stop codons by expanding aligned region in genome.
    
    Expand & Check strategy with STRICT FRAME verification:
    1. Extract genomic region with padding (max 20% of ref length)
    2. Orient sequence (reverse complement if negative strand)
    3. Search for start codon upstream (if check_start=True) IN FRAME (steps of 3)
    4. Search for stop codon downstream (if check_stop=True) IN FRAME (steps of 3)
    
    Args:
        genome_fasta_path: Path to genome FASTA file
        chromosome: Chromosome/contig ID
        start: Alignment start (1-based)
        end: Alignment end (1-based)
        strand: '+' or '-'
        reference_length_aa: Length of reference protein (for window calculation)
        check_start: Search for start codon (default: True)
        check_stop: Search for stop codon (default: True)
    
    Returns:
        Dict with 'has_start', 'has_stop', 'start_codon', 'stop_codon'
    """
    # Load genome sequence
    genome_seqs = {rec.id: rec.seq for rec in SeqIO.parse(genome_fasta_path, "fasta")}
    
    if chromosome not in genome_seqs:
        logger.warning(f"Chromosome {chromosome} not found in genome")
        return {"has_start": False, "has_stop": False, "start_codon": None, "stop_codon": None}
    
    genome_seq = genome_seqs[chromosome]
    
    # Calculate dynamic padding (max 20% of reference length)
    # reference_length_aa * 3 = reference_length_nt
    # 20% of reference_length_nt
    max_padding_nt = int((reference_length_aa * 3) * 0.20)
    padding = max_padding_nt
    
    # Convert to 0-based
    start_0 = start - 1
    end_0 = end
    
    # Adjust coordinates with padding
    search_start = max(0, start_0 - padding)
    search_end = min(len(genome_seq), end_0 + padding)
    
    # Extract extended region
    raw_seq = genome_seq[search_start:search_end]
    
    # Orient sequence (reverse complement if negative strand)
    if strand == '-':
        raw_seq = raw_seq.reverse_complement()
    
    dna_str = str(raw_seq).upper()
    
    # Calculate relative alignment position within extracted string
    if strand == '+':
        rel_align_start = start_0 - search_start
        rel_align_end = end_0 - search_start
    else:
        # On reverse strand, after reverse_complement: original 'end' becomes 5' start
        rel_align_start = search_end - end_0
        rel_align_end = search_end - start_0
    
    # Valid bacterial start codons (E. coli)
    valid_starts = ["ATG", "GTG", "TTG"]
    stop_codons = ["TAA", "TAG", "TGA"]
    
    # --- START CODON SEARCH (UPSTREAM) ---
    has_start = False
    start_codon = None
    
    if check_start:
        upstream_region = dna_str[0:rel_align_start]
        
        # Iterate backwards searching for start codon IN FRAME (steps of 3)
        # We start from the codon immediately preceding the alignment start
        # i.e., rel_align_start - 3, rel_align_start - 6, etc.
        for i in range(len(upstream_region) - 3, -1, -3):
            codon = upstream_region[i:i+3]
            if len(codon) < 3:
                continue
            
            if codon in valid_starts:
                has_start = True
                start_codon = codon
                break
            
            if codon in stop_codons:
                # Stop before start = truncated gene / pseudogene
                break
    
    # --- STOP CODON SEARCH (DOWNSTREAM) ---
    has_stop = False
    stop_codon = None
    
    if check_stop:
        downstream_region = dna_str[rel_align_end:]
        
        # Iterate forward searching for stop codon IN FRAME (steps of 3)
        # We start from the codon immediately following the alignment end
        # i.e., 0, 3, 6... relative to downstream_region start
        for i in range(0, len(downstream_region) - 2, 3):
            codon = downstream_region[i:i+3]
            if len(codon) < 3:
                continue
            
            if codon in stop_codons:
                has_stop = True
                stop_codon = codon
                break
            
            # Note: We don't break on start codons downstream, as they might be Methionines inside the tail
    
    return {
        "has_start": has_start,
        "has_stop": has_stop,
        "start_codon": start_codon,
        "stop_codon": stop_codon
    }


def analyze_protein_for_disablements(
    protein_hit_group: ProteinHitGroup
) -> DisablementCounts:
    """
    Analyze HSP group for inactivating mutations.
    
    Args:
        protein_hit_group: Protein HSP group
    
    Returns:
        Counters for different mutation types
    """
    logger.debug(f"Analyzing protein {protein_hit_group.protein_id}")
    
    counts = DisablementCounts()
    
    hsps = protein_hit_group.hsps
    
    if not hsps:
        return counts
    
    # Concatenate all HSP sequences for global analysis
    all_query_seq = ""
    all_subject_seq = ""
    
    for hsp in hsps:
        all_query_seq += hsp.qseq
        all_subject_seq += hsp.sseq
    
    # 1. Non-synonymous substitutions
    counts.non_synonymous_substitutions = count_non_synonymous_substitutions(
        all_query_seq, all_subject_seq
    )
    
    # 2. In-frame indels
    counts.in_frame_indels = count_in_frame_indels(
        all_query_seq, all_subject_seq
    )
    
    # 3. Frameshifts
    counts.frameshifts = count_frameshifts(hsps)
    
    # 4. Missing start codon (will be verified genomically in annotate_pseudogenes)
    counts.missing_start_codon = 1 if check_missing_start_codon(all_subject_seq) else 0
    
    # 5. Missing stop codon (will be verified genomically in annotate_pseudogenes)
    counts.missing_stop_codon = 1 if check_missing_stop_codon(all_subject_seq) else 0
    
    # 6. Premature stop codons
    counts.premature_stop_codons = count_premature_stop_codons(all_subject_seq)
    
    logger.debug(f"{protein_hit_group.protein_id}: {counts.total_disablements} mutations")
    
    return counts


def classify_as_pseudogene(
    disablements: DisablementCounts,
    min_disablements: int = 1
) -> bool:
    """
    Classifica se uma região é um pseudogene baseado no número de mutações.
    
    CRITÉRIOS ATUALIZADOS:
    - Stop codons prematuros são forte evidência de pseudogene
    - Frameshifts são forte evidência de pseudogene
    - Ausência de start codon (verificada genomicamente) é forte evidência
    - Ausência de stop codon (verificada genomicamente) é forte evidência
    - Tamanho <80% ou >120% da proteína de referência é forte evidência (Regra dos 20%)
    - Substituições e indels sozinhos NÃO indicam pseudogene (podem ser variação normal)
    
    Args:
        disablements: Contadores de mutações
        min_disablements: Número mínimo de mutações para classificar como pseudogene
    
    Returns:
        True se for classificado como pseudogene, False caso contrário
    """
    # Evidências fortes de pseudogene
    strong_evidence = (
        disablements.premature_stop_codons +
        disablements.frameshifts +
        disablements.missing_start_codon +
        disablements.missing_stop_codon +
        disablements.size_mismatch
    )
    
    # Classifica como pseudogene apenas se houver evidência forte
    is_pseudogene = strong_evidence >= min_disablements
    
    logger.debug(
        f"Classificação: strong_evidence={strong_evidence}, "
        f"is_pseudogene={is_pseudogene}"
    )
    
    return is_pseudogene


def annotate_pseudogenes(
    region_annotations: List[RegionAnnotation],
    genome_fasta_path: Path,
    output_path: Path,
    min_disablements: int = 1,
    padding: int = 150
) -> List[PseudogeneAnnotation]:
    """
    Anota regiões genômicas como potenciais pseudogenes.
    
    Usa verificação genômica "Expand & Check" para detectar ausência
    de start/stop codons na sequência de DNA real.
    
    Args:
        region_annotations: Lista de anotações de região com melhores proteínas
        genome_fasta_path: Caminho para o arquivo FASTA do genoma
        min_disablements: Número mínimo de mutações para classificar como pseudogene
        output_path: Caminho para salvar a lista de anotações de pseudogenes 
        padding: Nucleotídeos extras para buscar start/stop (padrão: 150bp)
    
    Returns:
        Lista de anotações de pseudogenes
    """
    logger.debug(f"Detecting pseudogenes in {len(region_annotations)} regions...")
    
    pseudogene_annotations = []
    
    for region_ann in region_annotations:
        # Analyze best protein for inactivating mutations
        disablements = analyze_protein_for_disablements(region_ann.best_protein)
        
        # --- SMART START/STOP VERIFICATION ---
        # Step 1: Quick check in alignment (if Methionine present, start codon exists)
        # Step 2: Genomic verification only if alignment check inconclusive
        
        hsps = region_ann.best_protein.hsps
        sorted_hsps = sorted(hsps, key=lambda h: h.qstart)
        first_hsp = sorted_hsps[0]
        last_hsp = sorted_hsps[-1]
        
        # Check start: if alignment has Methionine at beginning (any qstart), start codon exists
        first_qseq_clean = first_hsp.qseq.lstrip('-')
        has_start_in_alignment = (
            len(first_qseq_clean) > 0 and
            first_qseq_clean[0] == 'M'
            # Note: qstart position doesn't matter - divergent N-terminals are normal
        )
        
        # Check stop: if alignment reaches end of query, stop codon exists
        has_stop_in_alignment = (last_hsp.qend == last_hsp.qlen)
        
        # Genomic verification (only if alignment check failed)
        need_genomic_check = not has_start_in_alignment or not has_stop_in_alignment
        
        if need_genomic_check:
            genomic_check = check_start_stop_codons_genomic(
                genome_fasta_path,
                region_ann.region.chromosome,
                region_ann.region.start,
                region_ann.region.end,
                region_ann.region.strand,
                reference_length_aa=region_ann.best_protein.hsps[0].qlen,
                check_start=not has_start_in_alignment,
                check_stop=not has_stop_in_alignment
            )
            
            # Update based on genomic check (only if alignment didn't confirm)
            if not has_start_in_alignment:
                disablements.missing_start_codon = 0 if genomic_check["has_start"] else 1
            else:
                disablements.missing_start_codon = 0  # Confirmed in alignment
            
            if not has_stop_in_alignment:
                disablements.missing_stop_codon = 0 if genomic_check["has_stop"] else 1
            else:
                disablements.missing_stop_codon = 0  # Confirmed in alignment
        else:
            # Both confirmed in alignment
            disablements.missing_start_codon = 0
            disablements.missing_stop_codon = 0
        
        # --- SIZE-BASED PSEUDOGENE DETECTION (20% RULE) ---
        # TEMPORARILY DISABLED - Under investigation
        # TODO: Review biological assumptions and implementation
        """
        # Apply only if both start and stop codons are present (per specification)
        has_both_codons = (disablements.missing_start_codon == 0 and 
                          disablements.missing_stop_codon == 0)
        
        if has_both_codons:
            # Calculate candidate length from query coverage (handles overlapping HSPs)
            hsps = region_ann.best_protein.hsps
            min_qstart = min(hsp.qstart for hsp in hsps)
            max_qend = max(hsp.qend for hsp in hsps)
            candidate_length_aa = max_qend - min_qstart + 1
            
            # Reference protein length (all HSPs share same qlen)
            reference_length_aa = hsps[0].qlen
            length_check = check_length_viability(candidate_length_aa, reference_length_aa)
            
            if length_check["status"] == "FAIL":
                disablements.size_mismatch = 1
                logger.debug(
                    f"Size mismatch: candidate={length_check['candidate_len']}aa, "
                    f"reference={length_check['ref_len']}aa, "
                    f"ratio={length_check['ratio']:.2f} - {length_check['classification']}"
                )
            else:
                disablements.size_mismatch = 0
                logger.debug(
                    f"Size check passed: candidate={length_check['candidate_len']}aa, "
                    f"reference={length_check['ref_len']}aa, ratio={length_check['ratio']:.2f}"
                )
        else:
            # Skip size check if start/stop codons are missing
            disablements.size_mismatch = 0
        """
        # Force size_mismatch to 0 while disabled
        disablements.size_mismatch = 0
        
        # Check if region is a Small ORF
        # Condition: Genomic length < 300 nt AND Alignment length < 100 aa
        # (If EITHER is large enough, it is NOT a small ORF)
        genomic_len = region_ann.region.length
        alignment_len = region_ann.best_protein.total_length
        
        is_normal_size = (genomic_len >= MIN_REGION_SIZE_NT) or (alignment_len >= MIN_ALIGNMENT_SIZE_AA)
        is_small_orf = not is_normal_size
        
        # Classify as pseudogene based on mutations and size ratio
        # Small ORFs are classified normally (likely pseudogenes if truncated, functional if full-length small proteins)
        is_pseudogene = classify_as_pseudogene(disablements, min_disablements)
        
        pseudogene_ann = PseudogeneAnnotation(
            region_annotation=region_ann,
            disablements=disablements,
            is_pseudogene=is_pseudogene,
            is_small_orf=is_small_orf
        )
        
        pseudogene_annotations.append(pseudogene_ann)
    
    num_pseudogenes = sum(1 for ann in pseudogene_annotations if ann.is_pseudogene)
    logger.debug(f"Pseudogenes: {num_pseudogenes}/{len(pseudogene_annotations)}")

    logger.debug(f"Saving {len(pseudogene_annotations)} pseudogene annotations to: {output_path}")

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open('w') as f:
        for ann in pseudogene_annotations:
            f.write(json.dumps(ann.to_json_dict()) + '\n')  # JSON Lines: 1 objeto por linha


    # data = [ann.to_dict() for ann in pseudogene_annotations]
    # df = pd.DataFrame(data)
    # output_path.parent.mkdir(parents=True, exist_ok=True)
    # df.to_csv(output_path, sep='\t', index=False)

    return pseudogene_annotations


def calculate_subject_coverage_nt(hsps: List[BlastHit]) -> int:
    """
    Calcula a cobertura no SUBJECT (Genoma) em nucleotídeos, excluindo gaps.
    
    Args:
        hsps: Lista de HSPs
        
    Returns:
        Total de nucleotídeos cobertos pelo alinhamento
    """
    if not hsps:
        return 0
    
    # Coletar intervalos no genoma (Subject)
    # Nota: tblastn pode ter sstart > send se for na fita reversa
    intervals = []
    for hsp in hsps:
        start = min(hsp.sstart, hsp.send)
        end = max(hsp.sstart, hsp.send)
        intervals.append((start, end))
    
    # Ordenar por posição inicial
    intervals.sort()
    
    # Fundir intervalos sobrepostos
    merged = []
    for start, end in intervals:
        if not merged or start > merged[-1][1] + 1:
            merged.append([start, end])
        else:
            merged[-1][1] = max(merged[-1][1], end)
    
    # Somar comprimentos
    total = sum(end - start + 1 for start, end in merged)
    
    return total


def load_region_annotations_from_json(annotations_path: Path) -> List[RegionAnnotation]:
    """
    Carrega RegionAnnotations de um arquivo JSON Lines.

    Args:
        annotations_path: Caminho para o arquivo .jsonl gerado por annotate_regions_with_best_proteins

    Returns:
        Lista de RegionAnnotation completamente reconstruída
    """
    annotations = []
    with annotations_path.open('r') as f:
        for line in f:
            line = line.strip()
            if line:
                d = json.loads(line)
                annotations.append(RegionAnnotation.from_dict(d))

    logger.debug(f"Loaded {len(annotations)} region annotations from {annotations_path}")
    return annotations

import sys
if __name__ == "__main__":
    region_annotations = load_region_annotations_from_json(Path(sys.argv[1]))
    annotate_pseudogenes(region_annotations, Path(sys.argv[2]), Path(sys.argv[4]), min_disablements=int(sys.argv[3]))