process BedtoolsMergeBlastHits {
    publishDir "${params.output_merge}", mode: "copy"

    input:
    path filtered_hits_path
    val output_dir
    val genome_id

    output:
    path "${output_dir}/blast_hits.bed"
    path "${output_dir}/merged_regions.bed"

    script:
    """
    python3 ${workflow.projectDir}/../src/bedtools_merge.py ${filtered_hits_path} ${output_dir} ${genome_id}
    """
}