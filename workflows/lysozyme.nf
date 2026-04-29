include { RunBlastPipeline } from "../modules/local/run_blast_pipeline.nf";
include { RunFilteringStep } from "../modules/local/run_filtering_step.nf";
include { RealignFilteredHitsParallel } from "../modules/local/realign_filtered_hits_parallel.nf";
include { FilterSSearchByEValue } from "../modules/local/filter_ssearch_by_evalue.nf";

workflow {
    main:
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
            log.error "No hits passed filtering"
            return false
        }
        else {
            return true
        }
    }

    realign_filtered_hits_out = RealignFilteredHitsParallel(
        filtering_step_out,
        params.lysozyme_fasta,
        params.genome_fasta,
        params.output_ssearch,
        params.ssearch_path
    )

    FilterSSearchByEValue(realign_filtered_hits_out)
}