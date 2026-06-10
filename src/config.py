"""
Configurações e constantes do pipeline de anotação de lisozimas.
Define todos os parâmetros utilizados no pipeline para evitar hardcoding.
"""

from dataclasses import dataclass
from typing import Dict, List


# ==================== CONSTANTES DE FILTRAGEM ====================

# Threshold de identidade percentual mínima para filtragem inicial
MIN_IDENTITY_THRESHOLD: float = 20.0

# Score mínimo BLOSUM62 para filtragem inicial
MIN_BLOSUM62_SCORE: int = 120

# E-value máximo para busca BLAST
MAX_EVALUE: float = 1.0

# E-value threshold for SSEARCH filtering (Filter 2)
SSEARCH_EVALUE_THRESHOLD: float = 1e-6  # 10e-7 from specification


# ==================== PARÂMETROS DO BLAST ====================

@dataclass
class BlastParameters:
    """Parâmetros para execução do tblastn."""

    word_size: int = 3
    gapopen: int = 11
    gapextend: int = 1
    matrix: str = "BLOSUM62"
    comp_based_stats: int = 2
    seg: str = "yes"
    soft_masking: bool = True
    lcase_masking: bool = True
    evalue: float = 1.0  # E-value threshold as per specification
    max_target_seqs: int = 1000000
    threshold: int = 13
    num_threads: int = 1  # Será configurado pelo módulo de dependências

    # Formato de saída customizado
    outfmt_fields: List[str] = None

    def __post_init__(self):
        if self.outfmt_fields is None:
            self.outfmt_fields = [
                "qseqid",    # Query sequence ID
                "qlen",      # Query sequence length
                "sseqid",    # Subject sequence ID
                "slen",      # Subject sequence length
                "qstart",    # Query alignment start
                "qend",      # Query alignment end
                "sstart",    # Subject alignment start
                "send",      # Subject alignment end
                "qseq",      # Query aligned sequence
                "sseq",      # Subject aligned sequence
                "evalue",    # E-value
                "bitscore",  # Bit score
                "score",     # Raw score
                "length",    # Alignment length
                "pident",    # Percentage of identical matches
                "nident",    # Number of identical matches
                "mismatch",  # Number of mismatches
                "positive",  # Number of positive-scoring matches
                "gapopen",   # Number of gap openings
                "gaps",      # Total number of gaps
                "ppos",      # Percentage of positive-scoring matches
                "sframe",    # Subject frame
                "sstrand",   # Subject strand
                "qcovs",     # Query coverage per subject
                "qcovhsp"    # Query coverage per HSP
            ]
    
    def get_outfmt_string(self) -> str:
        """Retorna a string de formato de saída para o comando BLAST."""
        fields_str = " ".join(self.outfmt_fields)
        return f"6 {fields_str}"
    
    def to_command_args(self) -> List[str]:
        """Converte os parâmetros em argumentos de linha de comando."""
        args = [
            "-word_size", str(self.word_size),
            "-gapopen", str(self.gapopen),
            "-gapextend", str(self.gapextend),
            "-matrix", self.matrix,
            "-comp_based_stats", str(self.comp_based_stats),
            "-seg", self.seg,
            "-evalue", str(self.evalue),
            "-max_target_seqs", str(self.max_target_seqs),
            "-threshold", str(self.threshold),
            "-num_threads", str(self.num_threads),
            "-outfmt", self.get_outfmt_string()
        ]
        
        if self.soft_masking:
            args.append("-soft_masking")
            args.append("true")
        
        if self.lcase_masking:
            args.append("-lcase_masking")
        
        return args


# ==================== PARÂMETROS DO SSEARCH ====================

@dataclass
class SSearchParameters:
    """Parâmetros para execução do ssearch36 (Smith-Waterman)."""
    
    best_hits: int = 1          # -b: número de melhores hits a mostrar
    best_alignments: int = 1    # -d: número de melhores alinhamentos
    shuffle_count: int = 500    # -k: número de permutações estatísticas
    gap_open_penalty: int = -11 # -f: penalidade de abertura de gap
    gap_extend_penalty: int = 1   # -g: penalidade de extensão de gap (positive value)
    output_format: str = "8C"   # -m: formato de saída (8C = tabular tipo BLAST com comentários)
    quiet_mode: bool = True     # -q: modo silencioso
    scoring_matrix: str = "BL62"  # -s: matriz de scoring
    statistic_calculation: int = 11  # -z: tipo de cálculo estatístico
    database_size: int = 10000  # -Z: tamanho aparente do banco de dados
    
    def to_command_args(self) -> List[str]:
        """Converte os parâmetros em argumentos de linha de comando."""
        args = [
            "-b", str(self.best_hits),
            "-d", str(self.best_alignments),
            "-k", str(self.shuffle_count),
            "-f", str(self.gap_open_penalty),
            "-g", str(self.gap_extend_penalty),
            "-m", self.output_format,
            "-s", self.scoring_matrix,
            "-z", str(self.statistic_calculation),
            "-Z", str(self.database_size)
        ]
        
        if self.quiet_mode:
            args.append("-q")
        
        return args


# ==================== PARÂMETROS DO BEDTOOLS ====================

@dataclass
class BedToolsParameters:
    """Parâmetros para execução do bedtools merge."""
    
    # Colunas para operações de agregação
    # -c: colunas a agregar (2=nome, 4=score, 5=strand)
    aggregate_columns: str = "2,4,4,4,5"
    
    # Operações de agregação para cada coluna
    # -o: count, mean, min, max, collapse
    aggregate_operations: str = "count,mean,min,max,collapse"
    
    def to_command_args(self) -> List[str]:
        """Converte os parâmetros em argumentos de linha de comando."""
        return [
            "-c", self.aggregate_columns,
            "-o", self.aggregate_operations
        ]


# ==================== CONSTANTES DE DETECÇÃO DE PSEUDOGENES ====================

# Caractere que representa gap em alinhamentos
GAP_CHAR: str = "-"

# Caractere que representa stop codon em sequências traduzidas
STOP_CODON_CHAR: str = "*"

# Aminoácido inicial esperado (Metionina)
START_CODON_AA: str = "M"

# Limites de tamanho para a Regra dos 20%
MIN_LENGTH_RATIO: float = 0.8  # 80% do tamanho da referência
MAX_LENGTH_RATIO: float = 1.2  # 120% do tamanho da referência

# Minimum region size for Small ORF classification (nucleotides)
MIN_REGION_SIZE_NT: int = 300  # Regions below this are classified as 'small ORF' (100 aa)
MIN_ALIGNMENT_SIZE_AA: int = 100  # Alignment length threshold for Small ORF check


# ==================== ESTRUTURAS DE DADOS ====================

@dataclass
class DisablementCounts:
    """Contadores de diferentes tipos de mutações que podem inativar um gene."""
    
    non_synonymous_substitutions: int = 0  # Substituições não-sinônimas
    in_frame_indels: int = 0               # Indels in-frame
    frameshifts: int = 0                   # Mudanças de frame
    missing_start_codon: int = 0           # Perda do start codon
    missing_stop_codon: int = 0            # Perda do stop codon
    premature_stop_codons: int = 0         # Stop codons prematuros
    size_mismatch: int = 0                 # Tamanho incompatível (Regra dos 20%)
    
    @property
    def total_disablements(self) -> int:
        """Retorna o total de mutações inativadoras."""
        return (
            self.non_synonymous_substitutions +
            self.in_frame_indels +
            self.frameshifts +
            self.missing_start_codon +
            self.missing_stop_codon +
            self.premature_stop_codons +
            self.size_mismatch
        )
    
    def to_dict(self) -> Dict[str, int]:
        """Converte os contadores em dicionário."""
        return {
            "non_synonymous_substitutions": self.non_synonymous_substitutions,
            "in_frame_indels": self.in_frame_indels,
            "frameshifts": self.frameshifts,
            "missing_start_codon": self.missing_start_codon,
            "missing_stop_codon": self.missing_stop_codon,
            "premature_stop_codons": self.premature_stop_codons,
            "size_mismatch": self.size_mismatch,
            "total_disablements": self.total_disablements
        }
    
    @classmethod
    def from_dict(cls, d: Dict) -> 'DisablementCounts':
        return cls(
            non_synonymous_substitutions=d['non_synonymous_substitutions'],
            in_frame_indels=d['in_frame_indels'],
            frameshifts=d['frameshifts'],
            missing_start_codon=d['missing_start_codon'],
            missing_stop_codon=d['missing_stop_codon'],
            premature_stop_codons=d['premature_stop_codons'],
            size_mismatch=d['size_mismatch'],
            # total_disablements ignorado — é @property calculada automaticamente
        )


# ==================== INSTÂNCIAS PADRÃO ====================

# Instâncias padrão dos parâmetros para uso no pipeline
DEFAULT_BLAST_PARAMS = BlastParameters()
DEFAULT_SSEARCH_PARAMS = SSearchParameters()
DEFAULT_BEDTOOLS_PARAMS = BedToolsParameters()
