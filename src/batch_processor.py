"""
Module for batch processing of multiple genomes.
Implements the controller that iterates over multiple FASTA files.
"""

import logging
import glob
import pandas as pd
from pathlib import Path
from typing import List, Tuple, Dict
from datetime import datetime

from src.blast_search import run_blast_pipeline_step
from src.blast_filter import run_filtering_step, parse_blast_output
from src.ssearch_realign import realign_filtered_hits_parallel
from src.bedtools_merge import merge_blast_hits, save_merged_regions
from src.score_density import annotate_regions_with_best_proteins, save_region_annotations
from src.pseudogene_detection import (
    annotate_pseudogenes, 
    save_pseudogene_annotations,
    save_coverage_statistics,
    generate_summary_report
)
from src.export_gff3 import export_to_gff3
from src.comparative_analysis import generate_comparative_report
from src.config import (
    DEFAULT_BLAST_PARAMS,
    DEFAULT_SSEARCH_PARAMS,
    DEFAULT_BEDTOOLS_PARAMS,
    MIN_IDENTITY_THRESHOLD,
    MIN_BLOSUM62_SCORE
)


logger = logging.getLogger(__name__)


class BatchProcessor:
    """Batch processor for multiple genomes."""
    
    def __init__(
        self,
        input_dir: Path,
        lysozymes_path: Path,
        output_dir: Path,
        dependencies: dict,
        min_identity: float = MIN_IDENTITY_THRESHOLD,
        min_score: int = MIN_BLOSUM62_SCORE,
        min_disablements: int = 1,
        min_coverage: float = 0.8,
        final_min_identity: float = 0.0,
        num_threads: int = None
    ):
        """
        Initialize batch processor.
        
        Args:
            input_dir: Directory with genome FASTA files
            lysozymes_path: FASTA file with reference lysozymes
            output_dir: Base output directory
            dependencies: Dictionary with tool paths
            min_identity: Minimum identity for filtering
            min_score: Minimum score for filtering
            min_disablements: Minimum mutations for pseudogene
            min_coverage: Minimum coverage ratio for reporting
            final_min_identity: Final minimum identity filter
            num_threads: Number of threads for parallelism
        """
        self.input_dir = Path(input_dir)
        self.lysozymes_path = Path(lysozymes_path)
        self.output_dir = Path(output_dir)
        self.dependencies = dependencies
        self.min_identity = min_identity
        self.min_score = min_score
        self.min_disablements = min_disablements
        self.min_coverage = min_coverage
        self.final_min_identity = final_min_identity
        self.num_threads = num_threads
        
        self.per_genome_dir = output_dir / "per_genome"
        self.comparative_dir = output_dir / "comparative_analysis"
        
    def discover_genomes(self) -> List[Tuple[str, Path]]:
        """
        Discover FASTA files in input directory.
        
        Returns:
            List of tuples (genome_id, genome_path)
        """
        logger.debug(f"Searching for genomes in: {self.input_dir}")
        
        # Search for .fasta, .fa files
        fasta_files = []
        for pattern in ['*.fasta', '*.fa', '*.fna']:
            fasta_files.extend(glob.glob(str(self.input_dir / pattern)))
        
        if not fasta_files:
            raise FileNotFoundError(f"Nenhum arquivo FASTA encontrado em: {self.input_dir}")
        
        genomes = []
        for fasta_path in fasta_files:
            path = Path(fasta_path)
            # Use filename (without extension) as genome_id
            genome_id = path.stem
            genomes.append((genome_id, path))
        
        logger.debug(f"{len(genomes)} genomes found")
        for genome_id, _ in genomes:
            logger.debug(f"  - {genome_id}")
        
        return genomes
    
    def process_single_genome(
        self,
        genome_id: str,
        genome_path: Path,
        genome_num: int = 0,
        total_genomes: int = 0
    ) -> pd.DataFrame:
        """
        Execute complete pipeline for a single genome.
        
        Args:
            genome_id: Genome identifier
            genome_path: Path to genome FASTA file
            genome_num: Current genome number (for progress display)
            total_genomes: Total number of genomes (for progress display)
        
        Returns:
            DataFrame with final results (pseudogene_annotations)
        """
        logger.info(f"[{genome_num}/{total_genomes}] {genome_id}...")
        
        # Create output directories for this genome
        genome_output_dir = self.per_genome_dir / genome_id
        blast_dir = genome_output_dir / "blast"
        ssearch_dir = genome_output_dir / "ssearch"
        merge_dir = genome_output_dir / "merge"
        final_dir = genome_output_dir / "final"
        
        for dir_path in [blast_dir, ssearch_dir, merge_dir, final_dir]:
            dir_path.mkdir(parents=True, exist_ok=True)
        
        try:
            # Configure number of threads for BLAST
            blast_threads = self.num_threads if self.num_threads is not None else self.dependencies['num_threads']
            DEFAULT_BLAST_PARAMS.num_threads = blast_threads
            
            # STEP 1: BLAST
            logger.debug(f"[{genome_id}] STEP 1: Homology Search (BLAST)")
            logger.debug(f"  Using {blast_threads} threads for BLAST")
            
            blast_output = run_blast_pipeline_step(
                genome_path,
                self.lysozymes_path,
                blast_dir,
                DEFAULT_BLAST_PARAMS,
                self.dependencies['makeblastdb'],
                self.dependencies['tblastn']
            )
            
            # Parse with genome_id
            blast_hits = parse_blast_output(blast_output, genome_id=genome_id)
            logger.debug(f"  BLAST hits: {len(blast_hits)}")
            
            # STEP 2: FILTERING
            logger.debug(f"[{genome_id}] STEP 2: Initial Filtering")
            filtered_output = blast_dir / "filtered_hits.tsv"
            filtered_hits = run_filtering_step(
                blast_output,
                filtered_output,
                self.min_identity,
                self.min_score,
                genome_id=genome_id
            )
            
            if not filtered_hits:
                logger.warning(f"[{genome_id}] No hits passed filtering")
                return pd.DataFrame()
            
            logger.debug(f"  Filtered hits: {len(filtered_hits)}")
            
            # STEP 3: SSEARCH (Parallel)
            logger.debug(f"[{genome_id}] STEP 3: SSEARCH Realignment (Parallel)")
            ssearch_alignments = realign_filtered_hits_parallel(
                filtered_hits=filtered_hits,
                query_fasta_path=self.lysozymes_path,
                genome_fasta_path=genome_path,
                output_dir=ssearch_dir,
                ssearch_path=self.dependencies['ssearch36'],
                num_threads=self.num_threads
            )
            
            # Filter by E-value
            evalue_threshold = 1e-7
            filtered_ssearch = {key: aln for key, aln in ssearch_alignments.items() 
                               if aln.evalue <= evalue_threshold}
            
            logger.debug(f"  SSEARCH realignments (E-value ≤ {evalue_threshold:.0e}): {len(filtered_ssearch)}")
            
            # Update hits based on SSEARCH
            # AND update their scores with SSEARCH scores
            ssearch_hit_keys = set(filtered_ssearch.keys())
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
            
            logger.debug(f"  Updated {len(filtered_hits_after_ssearch)} hits with SSEARCH scores")
            filtered_hits = filtered_hits_after_ssearch
            
            if not filtered_hits:
                logger.warning(f"[{genome_id}] No hits passed SSEARCH filtering")
                return pd.DataFrame()
            
            # STEP 4: SEGMENT MERGING
            logger.debug(f"[{genome_id}] STEP 4: Segment Merging (BEDTools)")
            merged_regions = merge_blast_hits(
                filtered_hits,
                merge_dir,
                DEFAULT_BEDTOOLS_PARAMS,
                self.dependencies['bedtools'],
                genome_id=genome_id
            )
            logger.debug(f"  Merged regions: {len(merged_regions)}")
            
            # STEP 5: SCORE DENSITY
            logger.debug(f"[{genome_id}] STEP 5: Score Density Calculation")
            region_annotations = annotate_regions_with_best_proteins(
                merged_regions,
                filtered_hits
            )
            logger.debug(f"  Annotated regions: {len(region_annotations)}")
            
            # STEP 6: PSEUDOGENE DETECTION
            logger.debug(f"[{genome_id}] STEP 6: Pseudogene Detection")
            pseudogene_annotations = annotate_pseudogenes(
                region_annotations,
                genome_path,
                self.min_disablements
            )
            
            # --- Coverage Filter ---
            if self.min_coverage > 0:
                logger.debug(f"[{genome_id}] Applying coverage filter: >= {self.min_coverage*100:.1f}%")
                original_count = len(pseudogene_annotations)
                
                filtered_annotations = []
                for ann in pseudogene_annotations:
                    hsps = ann.region_annotation.best_protein.hsps
                    if hsps:
                        min_qstart = min(hsp.qstart for hsp in hsps)
                        max_qend = max(hsp.qend for hsp in hsps)
                        coverage_len = max_qend - min_qstart + 1
                        ref_len = hsps[0].qlen
                        coverage_ratio = coverage_len / ref_len if ref_len > 0 else 0
                        
                        if coverage_ratio >= self.min_coverage:
                            filtered_annotations.append(ann)
                
                pseudogene_annotations = filtered_annotations
                logger.debug(f"  Filtered {original_count - len(pseudogene_annotations)} regions. Remaining: {len(pseudogene_annotations)}")

            # --- Final Identity Filter ---
            if self.final_min_identity > 0:
                # Convert fraction to percentage if necessary (e.g. 0.7 -> 70.0)
                threshold_pct = self.final_min_identity * 100 if self.final_min_identity <= 1.0 else self.final_min_identity
                
                logger.debug(f"[{genome_id}] Applying final identity filter: >= {threshold_pct:.1f}%")
                original_count = len(pseudogene_annotations)
                
                pseudogene_annotations = [
                    ann for ann in pseudogene_annotations 
                    if max(hsp.pident for hsp in ann.region_annotation.best_protein.hsps) >= threshold_pct
                ]
                logger.debug(f"  Filtered {original_count - len(pseudogene_annotations)} regions. Remaining: {len(pseudogene_annotations)}")
            
            num_pseudogenes = sum(1 for a in pseudogene_annotations if a.is_pseudogene)
            functional = len(pseudogene_annotations) - num_pseudogenes
            logger.debug(f"  Pseudogenes: {num_pseudogenes}")
            logger.debug(f"  Functional genes: {functional}")
            
            # Save individual results
            pseudogenes_output = final_dir / "pseudogene_annotations.tsv"
            save_pseudogene_annotations(pseudogene_annotations, pseudogenes_output)
            
            # Save coverage statistics
            coverage_stats_output = final_dir / "coverage_statistics.tsv"
            save_coverage_statistics(pseudogene_annotations, coverage_stats_output)

            # Generate and save summary report
            summary_report = generate_summary_report(pseudogene_annotations, self.min_coverage)
            summary_output = final_dir / "lysozyme_summary.txt"
            with open(summary_output, 'w') as f:
                f.write(summary_report)
            
            # Export GFF3
            gff3_output = final_dir / f"{genome_id}_annotations.gff3"
            export_to_gff3(pseudogene_annotations, genome_id, gff3_output)
            
            num_pseudogenes = sum(1 for a in pseudogene_annotations if a.is_pseudogene)
            
            # Convert to DataFrame
            df = pd.DataFrame([{
                'genome_id': genome_id,
                'chromosome': a.region_annotation.region.chromosome,
                'start': a.region_annotation.region.start,
                'end': a.region_annotation.region.end,
                'best_protein_id': a.region_annotation.best_protein.protein_id,
                'score_density': a.region_annotation.best_protein.score_density,
                'is_pseudogene': a.is_pseudogene,
                'is_small_orf': a.is_small_orf,
                'non_synonymous_substitutions': a.disablements.non_synonymous_substitutions,
                'in_frame_indels': a.disablements.in_frame_indels,
                'frameshifts': a.disablements.frameshifts,
                'premature_stop_codons': a.disablements.premature_stop_codons,
                'missing_start_codon': a.disablements.missing_start_codon,
                'missing_stop_codon': a.disablements.missing_stop_codon,
                'total_disablements': a.disablements.total_disablements
            } for a in pseudogene_annotations])
            
            logger.info(f"  ✓ {genome_id}: {len(pseudogene_annotations)} regions, {num_pseudogenes} pseudogenes")
            
            return df
            
        except Exception as e:
            logger.error(f"[{genome_id}] Error during processing: {e}", exc_info=True)
            # Save error to file
            error_log = genome_output_dir / "error.log"
            with open(error_log, 'w') as f:
                f.write(f"Error processing {genome_id}:\n{str(e)}\n")
            return pd.DataFrame()
    
    def run_batch(self) -> pd.DataFrame:
        """
        Process all genomes in batch.
        
        Returns:
            Aggregated DataFrame with all results
        """
        start_time = datetime.now()
        
        logger.debug("Starting batch processing")
        
        # Discover genomes
        genomes = self.discover_genomes()
        total_genomes = len(genomes)
        
        # Process each genome
        all_results = []
        successful = 0
        failed = 0
        
        for idx, (genome_id, genome_path) in enumerate(genomes, start=1):
            df_result = self.process_single_genome(
                genome_id,
                genome_path,
                genome_num=idx,
                total_genomes=total_genomes
            )
            
            if not df_result.empty:
                all_results.append(df_result)
                successful += 1
            else:
                failed += 1
        
        # Generate comparative report
        if all_results:
            logger.debug("Generating comparative report")
            
            self.comparative_dir.mkdir(parents=True, exist_ok=True)
            
            from src.comparative_analysis import aggregate_all_genomes
            aggregated_df = aggregate_all_genomes(all_results)
            generate_comparative_report(aggregated_df, self.comparative_dir)
        else:
            logger.warning("No results to generate comparative report")
            aggregated_df = pd.DataFrame()
        
        # Final summary
        end_time = datetime.now()
        elapsed = end_time - start_time
        
        logger.info(f"\nComplete: {len(genomes)} genomes processed, {successful} successful, {failed} failed")
        logger.debug(f"Total execution time: {elapsed}")
        logger.debug(f"Results saved in: {self.output_dir}")
        
        return aggregated_df
