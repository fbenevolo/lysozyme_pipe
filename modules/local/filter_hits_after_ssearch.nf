process FilterHitsAfterSSearch {
    publishDir "${params.output_ssearch}", mode: "copy"

    input:
    path filtered_hits_path
    path filtered_ssearch_path

    output:
    path "filtered_hits_after_ssearch.tsv"

    script:
    """
    python3 ${workflow.projectDir}/../src/filter_hits_after_ssearch.py ${filtered_hits_path} ${filtered_ssearch_path} filtered_hits_after_ssearch.tsv
    """
}