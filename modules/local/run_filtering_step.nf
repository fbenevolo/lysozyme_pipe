process RunFilteringStep {
    publishDir "${params.output_blast}", mode: "copy"
    
    input:
    path blast_output
    path filtered_hits_output

    output:
    path "filtered_hits.tsv"

    script:
    """
    python3 ${workflow.projectDir}/../src/blast_filter.py ${blast_output} ${filtered_hits_output}
    """
}