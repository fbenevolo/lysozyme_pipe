include { RunBlastPipeline } from "../modules/local/run_blast_pipeline.nf";
include { RunFilteringStep } from "../modules/local/run_filtering_step.nf";
include { RealignFilteredHitsParallel } from "../modules/local/realign_filtered_hits_parallel.nf";
include { FilterSSearchByEValue } from "../modules/local/filter_ssearch_by_evalue.nf";
include { FilterHitsAfterSSearch } from "../modules/local/filter_hits_after_ssearch.nf";
include { BedtoolsMergeBlastHits } from "../modules/local/bedtools_merge_blast_hits.nf"; 
include { BedtoolsSaveMergedRegions } from "../modules/local/bedtools_save_merged_regions.nf";
include { AnnotateRegionsWithBestProteins } from "../modules/local/annotate_regions_with_best_proteins.nf";
include { AnnotatePseudogenes } from "../modules/local/annotate_pseudogenes.nf";
include { ApplyCoverage } from "../modules/local/apply_coverage.nf";
include { ApplyFinalIdentity } from "../modules/local/apply_final_identity.nf";
include { SavePseudogeneAnnotations } from "../modules/local/save_pseudogene_annotations.nf";
include { ExportToGFF3 } from "../modules/local/export_to_gff3.nf";
include { SaveCoverageStatistics } from "../modules/local/save_coverage_statistics.nf";
include { GenerateSummaryReport } from "../modules/local/generate_summary_report.nf";
include { BatchProcessor } from "../modules/local/batch_processor.nf";

workflow {
    main:

    if (params.input_dir) {
        BatchProcessor(
            params.input_dir,
            params.lysozyme_fasta,
            params.output_dir,
            params.min_identity,
            params.min_score,
            params.min_disablements,
            params.min_coverage,
            params.final_min_identity,
            params.num_threads
        )
        return;
    }


    def genome_id = file(params.genome_fasta).baseName;

    
    blast_out = RunBlastPipeline(
        file(params.genome_fasta), 
        file(params.lysozyme_fasta)
        )
    
    log.info "[1/3] Running Blast Search..."

    filtering_step_out = RunFilteringStep(
        blast_out[0],
        "${params.output_blast}/filtered_hits.tsv",
        genome_id
    )

    filtering_step_out.filter { file -> 
        if (file.readLines().size() <= 1) {
            error("No hits passed filtering")
        }
        else {
            log.info "Filtering step completed"
        }
    }

    realign_filtered_hits_out = RealignFilteredHitsParallel(
        filtering_step_out,
        params.lysozyme_fasta,
        params.genome_fasta,
        params.output_ssearch,
        params.ssearch_path
    )

    filtered_ssearch_by_evalue_out = FilterSSearchByEValue(realign_filtered_hits_out)
    
    filtered_ssearch_by_evalue_out.subscribe { log.info "SSearch realignments completed" }

    filtered_hits_after_ssearch_out = FilterHitsAfterSSearch(
        filtering_step_out,
        filtered_ssearch_by_evalue_out,
    )
    filtered_hits_after_ssearch_out.filter { file -> 
        if (file.readLines().size() <= 1) {
            error("No hits passed SSEARCH filtering")
        }
        else {
            log.info "Updated ${file.readLines().size()} hits with SSEARCH scores"
            log.info "[2/3] Merging ${file.readLines().size()} regions"
        }
    }

    merge_blast_hits_out = BedtoolsMergeBlastHits(
        filtered_hits_after_ssearch_out, 
        "./",
        genome_id
    )

    merged_regions_out = BedtoolsSaveMergedRegions(
        merge_blast_hits_out[1],
        "./merged_regions.tsv",
        genome_id
    )

    annotate_regions_with_best_proteins_out = AnnotateRegionsWithBestProteins(
        merge_blast_hits_out[1],
        filtered_hits_after_ssearch_out,
        genome_id
    )

    pseudogene_annotation_out = AnnotatePseudogenes(
        annotate_regions_with_best_proteins_out,
        params.genome_fasta,
        params.min_disablements
    )
    
    pseudogene_annotation_with_coverage_out = ApplyCoverage(
        pseudogene_annotation_out,
        params.min_coverage
    )

    pseudogene_annotation_parcial = (params.min_coverage > 0) ? pseudogene_annotation_with_coverage_out : pseudogene_annotation_out

    pseudogene_annotation_with_final_identity = ApplyFinalIdentity(
        pseudogene_annotation_parcial,
        params.final_min_identity
    )

    pseudogene_annotation_final = (params.final_min_identity > 0) ? pseudogene_annotation_with_final_identity : pseudogene_annotation_parcial

    SavePseudogeneAnnotations(pseudogene_annotation_final)
    ExportToGFF3(pseudogene_annotation_final, genome_id)
    SaveCoverageStatistics(pseudogene_annotation_final)
    GenerateSummaryReport(pseudogene_annotation_final, params.min_coverage)
}