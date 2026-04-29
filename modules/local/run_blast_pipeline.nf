process RunBlastPipeline {
    publishDir "${params.output_blast}", mode: "copy"

    input:
    path genome_fasta
    path lysozyme_fasta

    output:
    path "blast_results.tsv"
    path "genome_db.*"

    script:
    """
    python3 ${workflow.projectDir}/../src/blast_search.py ${genome_fasta} ${lysozyme_fasta} .
    """
}