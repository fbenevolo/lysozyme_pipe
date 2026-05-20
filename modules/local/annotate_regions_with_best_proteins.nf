process AnnotateRegionsWithBestProteins {
    publishDir "${params.output_final_dir}", mode: "copy"

    input:
    path genomic_regions_path
    path filtered_hits_path
    val genome_id

    output:
    path "region_annotations.jsonl"

    script:
    """
    python3 -u ${workflow.projectDir}/../src/score_density.py ${genomic_regions_path} ${filtered_hits_path} ./region_annotations.jsonl ${genome_id}
    """
}