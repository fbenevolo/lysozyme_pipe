include { RunBlastPipeline } from "../modules/local/run_blast_pipeline.nf";
include { RunFilteringStep } from "../modules/local/run_filtering_step.nf";
include { RealignFilteredHitsParallel } from "../modules/local/realign_filtered_hits_parallel.nf";
include { FilterSSearchByEValue } from "../modules/local/filter_ssearch_by_evalue.nf";
include { FilterHitsAfterSSearch } from "../modules/local/filter_hits_after_ssearch.nf";
include { MergeBlastHits } from "../modules/local/merge_blast_hits.nf"; 

workflow {
    main:
    
    log.info "[1/3] Running Blast Search..."
    blast_out = RunBlastPipeline(
        file(params.genome_fasta), 
        file(params.lysozyme_fasta)
        )

    filtering_step_out = RunFilteringStep(
        blast_out[0],
        "${params.output_blast}/filtered_hits.tsv"
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

    MergeBlastHits(
        filtered_hits_after_ssearch_out, 
        "./merge/"
    )
}