process FilterSSearchByEValue {
    publishDir "${params.output_ssearch}", mode: "copy"

    input:
    path realigned_filtered_hits

    output:
    path "ssearch_realignments_filtered_by_evalue.tsv"

    script:
    """
    python3 ${workflow.projectDir}/../src/filter_ssearch.py ${realigned_filtered_hits} ssearch_realignments_filtered_by_evalue.tsv
    """

}