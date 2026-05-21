#!/usr/bin/env python3
"""
Main lysozyme annotation and pseudogene detection pipeline.
Orchestrates all pipeline steps in an integrated manner.
"""

import argparse
import logging
import sys
from pathlib import Path
from datetime import datetime

# Importa todos os módulos do pipeline
from src.config import (
    DEFAULT_BLAST_PARAMS,
    DEFAULT_SSEARCH_PARAMS,
    DEFAULT_BEDTOOLS_PARAMS,
    MIN_IDENTITY_THRESHOLD,
    MIN_BLOSUM62_SCORE
)
from src.dependencies import verify_and_install_dependencies
from src.blast_search import run_blast_pipeline_step
from src.blast_filter import run_filtering_step
from src.ssearch_realign import realign_filtered_hits
from src.bedtools_merge import merge_blast_hits
from src.score_density import annotate_regions_with_best_proteins
from src.pseudogene_detection import (
    annotate_pseudogenes,
    save_coverage_statistics,
    generate_summary_report
)
from save_pseudogene_annotations import save_pseudogene_annotations
from apply_coverage_filter import apply_coverage
from apply_final_identity_filter import apply_final_identity
from src.export_gff3 import export_to_gff3
from src.filter_ssearch import filter_ssearch_by_evalue
from src.filter_hits_after_ssearch import filter_hits_after_ssearch
from src.bedtools_save_merged_regions import save_merged_regions


def setup_logging(log_file: Path = None, verbose: bool = False) -> None:
    """
    Configure logging system.
    
    Args:
        log_file: File to save logs (optional)
        verbose: If True, displays DEBUG messages
    """
    log_level = logging.DEBUG if verbose else logging.INFO
    
    # Message format
    log_format = '%(asctime)s - %(name)s - %(levelname)s - %(message)s'
    date_format = '%Y-%m-%d %H:%M:%S'
    
    # Basic configuration
    handlers = [logging.StreamHandler(sys.stdout)]
    
    if log_file:
        log_file.parent.mkdir(parents=True, exist_ok=True)
        handlers.append(logging.FileHandler(log_file))
    
    logging.basicConfig(
        level=log_level,
        format=log_format,
        datefmt=date_format,
        handlers=handlers
    )


def parse_arguments() -> argparse.Namespace:
    """
    Parse command line arguments.
    
    Returns:
        Namespace with parsed arguments
    """
    parser = argparse.ArgumentParser(
        description='Pipeline de anotação de lisozimas e detecção de pseudogenes em genomas de E. coli',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Exemplos de uso:
  
  # Execução básica
  python pipeline.py -g genome.fasta -l lysozymes.fasta -o results/
  
  # Com parâmetros customizados
  python pipeline.py -g genome.fasta -l lysozymes.fasta -o results/ \\
                     --min-identity 25 --min-score 150 -v
  
  # Apenas etapas específicas
  python pipeline.py -g genome.fasta -l lysozymes.fasta -o results/ \\
                     --skip-ssearch
        """
    )
    
    # Argumentos obrigatórios
    required = parser.add_argument_group('argumentos obrigatórios')
    
    # Modo single genome OU modo lote
    input_group = parser.add_mutually_exclusive_group(required=True)
    input_group.add_argument(
        '-g', '--genome',
        type=Path,
        metavar='FASTA',
        help='FASTA file of a single E. coli genome'
    )
    input_group.add_argument(
        '--input-dir',
        type=Path,
        metavar='DIR',
        help='Directory with multiple FASTA files (batch mode)'
    )
    
    required.add_argument(
        '-l', '--lysozymes',
        type=Path,
        required=True,
        metavar='FASTA',
        help='FASTA file with reference lysozymes (UniProtKB/Swiss-Prot)'
    )
    required.add_argument(
        '-o', '--output',
        type=Path,
        required=True,
        metavar='DIR',
        help='Output directory for results'
    )
    
    # Filtering parameters
    filtering = parser.add_argument_group('filtering parameters')
    filtering.add_argument(
        '--min-identity',
        type=float,
        default=MIN_IDENTITY_THRESHOLD,
        metavar='FLOAT',
        help=f'Minimum percent identity (default: {MIN_IDENTITY_THRESHOLD})'
    )
    filtering.add_argument(
        '--min-score',
        type=int,
        default=MIN_BLOSUM62_SCORE,
        metavar='INT',
        help=f'Minimum BLOSUM62 score (default: {MIN_BLOSUM62_SCORE})'
    )
    filtering.add_argument(
        '--min-disablements',
        type=int,
        default=1,
        metavar='INT',
        help='Minimum number of mutations to classify as pseudogene (default: 1)'
    )
    filtering.add_argument(
        '--min-coverage',
        type=float,
        default=0.8,
        metavar='FLOAT',
        help='Minimum coverage ratio for reporting (default: 0.8)'
    )
    filtering.add_argument(
        '--final-min-identity',
        type=float,
        default=0.0,
        metavar='FLOAT',
        help='Final minimum identity filter (default: 0.0 - disabled)'
    )
    
    # Execution options
    execution = parser.add_argument_group('execution options')
    execution.add_argument(
        '--num-threads',
        type=int,
        default=None,
        metavar='INT',
        help='Número de threads para SSEARCH paralelo (padrão: CPU_COUNT-1)'
    )
    execution.add_argument(
        '-v', '--verbose',
        action='store_true',
        help='Exibe mensagens de debug detalhadas'
    )
    execution.add_argument(
        '--log-file',
        type=Path,
        metavar='FILE',
        help='Arquivo para salvar logs (opcional)'
    )
    
    return parser.parse_args()


def validate_inputs(genome_path: Path, lysozymes_path: Path) -> None:
    """
    Valida os arquivos de entrada.
    
    Args:
        genome_path: Caminho para o arquivo do genoma
        lysozymes_path: Caminho para o arquivo de lisozimas
    
    Raises:
        FileNotFoundError: Se algum arquivo não existir
    """
    if not genome_path.exists():
        raise FileNotFoundError(f"Arquivo de genoma não encontrado: {genome_path}")
    
    if not lysozymes_path.exists():
        raise FileNotFoundError(f"Arquivo de lisozimas não encontrado: {lysozymes_path}")


def run_pipeline(
    genome_fasta: Path,
    lysozyme_fasta: Path,
    output_dir: Path,
    genome_id: str = "unknown",
    min_identity: float = MIN_IDENTITY_THRESHOLD,
    min_score: int = MIN_BLOSUM62_SCORE,
    min_disablements: int = 1,
    min_coverage: float = 0.8,
    final_min_identity: float = 0.0,
    num_threads: int = None,
    deps: dict = None
) -> None:
    """
    Executa o pipeline completo de anotação de lisozimas.
    
    Args:
        genome_fasta: Arquivo FASTA do genoma
        lysozyme_fasta: Arquivo FASTA com lisozimas de referência
        output_dir: Output directory
        genome_id: Genome identifier
        min_identity: Minimum identity for filtering
        min_score: Minimum score for filtering
        min_disablements: Minimum mutations to classify as pseudogene
        min_coverage: Minimum coverage ratio for reporting
        final_min_identity: Final minimum identity filter
        num_threads: Number of threads for parallelization
        deps: Dependencies dictionary from verify_and_install_dependencies()
    """
    logger = logging.getLogger(__name__)
    
    # If deps not provided, verify dependencies (for standalone usage)
    if deps is None:
        deps = verify_and_install_dependencies()
    
    # Configure BLAST threads
    from src.config import DEFAULT_BLAST_PARAMS
    blast_threads = num_threads if num_threads is not None else deps['num_threads']
    DEFAULT_BLAST_PARAMS.num_threads = blast_threads
    # Create output directories
    output_dir.mkdir(parents=True, exist_ok=True)
    blast_dir = output_dir / "blast"
    merge_dir = output_dir / "merge"
    final_dir = output_dir / "final"
    
    # STEP 1: BLAST search
    logger.info(f"[1/3] BLAST search ({blast_threads} threads)...")
    
    blast_output = run_blast_pipeline_step(
        genome_fasta,
        lysozyme_fasta,
        blast_dir,
        DEFAULT_BLAST_PARAMS,
        deps['makeblastdb'],
        deps['tblastn']
    )
    
    # STEP 2: Quality filtering (silent)
    filtered_output = blast_dir / "filtered_hits.tsv"
    filtered_hits = run_filtering_step(
        blast_output,
        filtered_output,
        min_identity,
        min_score,
        genome_id=genome_id
    )

    if not filtered_hits:
        logger.error("No hits passed filtering")
        return

    
    # STEP 3: SSEARCH realignment (silent)
    ssearch_dir = output_dir / "ssearch"
    # ssearch_dir.mkdir(parents=True, exist_ok=True)

    from src.ssearch_realign import realign_filtered_hits_parallel

    ssearch_alignments = realign_filtered_hits_parallel(
        filtered_hits=filtered_hits,
        query_fasta_path=lysozyme_fasta,
        genome_fasta_path=genome_fasta,
        output_dir=ssearch_dir,
        ssearch_path=deps['ssearch36'],
        num_threads=num_threads
    )

    # Filter realignments by E-value
    ssearch_output = ssearch_dir / "ssearch_realignments_filtered_by_evalue.tsv"
    filtered_ssearch = filter_ssearch_by_evalue(ssearch_alignments, ssearch_output)

    logger.debug(f"SSEARCH realignments complete: {len(ssearch_alignments)}")
    
    logger.debug(f"  - Before filtering: {len(ssearch_alignments)} realignments")
    logger.debug(f"  - After filtering: {len(filtered_ssearch)} realignments")
    
    ssearch_hit_keys = set(filtered_ssearch.keys())
    filtered_hits_after_ssearch = []

    # Update filtered hits to use only those that passed SSEARCH
    # AND update their scores with SSEARCH scores
    filtered_hits_after_ssearch_output = ssearch_dir / "filtered_hits_after_ssearch.tsv"
    filtered_hits_after_ssearch = filter_hits_after_ssearch(filtered_hits, filtered_ssearch, ssearch_hit_keys, filtered_hits_after_ssearch_output)
    if not filtered_hits_after_ssearch:
        logger.error("No hits passed SSEARCH filtering")
        return
    
    logger.debug(f"Updated {len(filtered_hits_after_ssearch)} hits with SSEARCH scores")
    filtered_hits = filtered_hits_after_ssearch
    
    # STEP 2: Region merging
    logger.info(f"[2/3] Merging {len(filtered_hits)} regions...")
    merged_regions = merge_blast_hits(
        filtered_hits,
        merge_dir,
        DEFAULT_BEDTOOLS_PARAMS,
        deps['bedtools'],
        genome_id=genome_id
    )
    
    
    merged_output = merge_dir / "merged_regions.tsv"
    save_merged_regions(merged_regions, merged_output)
    

    # Score density calculation (silent)
    annotations_output = final_dir / "region_annotations.jsonl"
    region_annotations = annotate_regions_with_best_proteins(
        merged_regions,
        filtered_hits,
        annotations_output
    )

    # STEP 3: Pseudogene detection
    logger.info(f"[3/3] Analyzing {len(region_annotations)} regions for pseudogenes...")
    
    initial_annotations_path = final_dir / "initial_pseudogene_annotations.jsonl"
    pseudogene_annotations = annotate_pseudogenes(
        region_annotations,
        genome_fasta,
        initial_annotations_path,
        min_disablements
    )
    
    # --- Coverage Filter ---
    if min_coverage > 0:
        logger.info(f"Applying coverage filter: >= {min_coverage*100:.1f}%")

        original_count = len(pseudogene_annotations)
        annotations_with_coverage_path = final_dir / "pseudogene_annotations_with_coverage.jsonl"
        filtered_annotations = apply_coverage(pseudogene_annotations, min_coverage, annotations_with_coverage_path)
        
        pseudogene_annotations = filtered_annotations
        logger.info(f"  Filtered {original_count - len(pseudogene_annotations)} regions. Remaining: {len(pseudogene_annotations)}")

    
    # --- Final Identity Filter ---
    if final_min_identity > 0:
        # Convert fraction to percentage if necessary (e.g. 0.7 -> 70.0)
        # If user provided > 1.0, assume it's already percentage
        threshold_pct = final_min_identity * 100 if final_min_identity <= 1.0 else final_min_identity
        
        logger.info(f"Applying final identity filter: >= {threshold_pct:.1f}%")
        original_count = len(pseudogene_annotations)

        annotations_with_final_identity_path = final_dir / "pseudogene_annotations_with_final_identity.jsonl"
        filtered_annotations = apply_final_identity(pseudogene_annotations, final_min_identity, annotations_with_final_identity_path)
        pseudogene_annotations = filtered_annotations
        
        logger.info(f"  Filtered {original_count - len(pseudogene_annotations)} regions. Remaining: {len(pseudogene_annotations)}")
    

    pseudogenes_output = final_dir / "pseudogene_annotations_final.tsv"
    logger.info(f"Saving {len(pseudogene_annotations)} pseudogene annotations to: {pseudogenes_output}")
    save_pseudogene_annotations(pseudogene_annotations, pseudogenes_output)
    logger.info("Pseudogene annotations saved successfully")
    
    
    # Export GFF3
    gff3_output = final_dir / "lysozyme_annotations.gff3"
    export_to_gff3(pseudogene_annotations, genome_id, gff3_output)
    
    '''
    # Save coverage statistics for detailed analysis
    coverage_output = final_dir / "coverage_statistics.tsv"
    save_coverage_statistics(pseudogene_annotations, coverage_output)
    
    # Generate summary report
    summary = generate_summary_report(pseudogene_annotations, min_coverage)
    print(summary)
    
    report_file = final_dir / "summary_report.txt"
    with open(report_file, 'w') as f:
        f.write(summary)
    
    num_pseudogenes = sum(1 for ann in pseudogene_annotations if ann.is_pseudogene)
    logger.info(f"Complete: {len(pseudogene_annotations)} regions, {num_pseudogenes} pseudogenes")
    '''

def main():
    """Função principal do pipeline."""
    # Parse dos argumentos
    args = parse_arguments()
    
    # Configura logging
    log_file = args.log_file if args.log_file else args.output / "pipeline.log"
    setup_logging(log_file, args.verbose)
    
    logger = logging.getLogger(__name__)
    
    try:
        # Verifica dependências
        deps = verify_and_install_dependencies()
        # Decide modo de execução: single genome vs. batch
        if args.input_dir:
            # ========== MODO LOTE ==========
            logger.info("Mode: BATCH (multiple genomes)")
            
            from src.batch_processor import BatchProcessor

            processor = BatchProcessor(
                input_dir=args.input_dir,
                lysozymes_path=args.lysozymes,
                output_dir=args.output,
                dependencies=deps,
                min_identity=args.min_identity,
                min_score=args.min_score,
                min_disablements=args.min_disablements,
                min_coverage=args.min_coverage,
                final_min_identity=args.final_min_identity,
                num_threads=args.num_threads
            )
            
            start_time = datetime.now()
            processor.run_batch()
            end_time = datetime.now()
        else:
            # ========== MODO SINGLE GENOME ==========
            logger.info("Mode: SINGLE GENOME")
            # Valida entradas
            validate_inputs(args.genome, args.lysozymes)
            
            # Extrai genome_id do nome do arquivo
            genome_id = args.genome.stem

            # Executa pipeline
            start_time = datetime.now()
            
            run_pipeline(
                genome_fasta=args.genome,
                lysozyme_fasta=args.lysozymes,
                output_dir=args.output,
                genome_id=genome_id,
                min_identity=args.min_identity,
                min_score=args.min_score,
                min_disablements=args.min_disablements,
                min_coverage=args.min_coverage,
                final_min_identity=args.final_min_identity,
                num_threads=args.num_threads,
                deps=deps  # Pass deps to avoid re-verification
            )
            end_time = datetime.now()
        
        elapsed = end_time - start_time
        
        logger.debug(f"Total execution time: {elapsed}")
        logger.debug("Pipeline completed successfully!")

    except Exception as e:
        logger.error(f"Erro durante execução do pipeline: {e}", exc_info=True)
        sys.exit(1)
        
    except Exception as e:
        logger.error(f"Erro durante execução do pipeline: {e}", exc_info=True)
        sys.exit(1)

if __name__ == "__main__":
    main()
