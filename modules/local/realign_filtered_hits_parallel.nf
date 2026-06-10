process RealignFilteredHitsParallel {
    publishDir "${params.output_ssearch}", mode: "copy"
    
    input:
    path filtered_hits
    path lysozyme_fasta
    path genome_fasta
    val output_dir
    val ssearch_path

    output:
    path "ssearch_realignments.tsv"

    script:
    """
    python3 ${workflow.projectDir}/../src/ssearch_realign.py ${filtered_hits} ${lysozyme_fasta} ${genome_fasta} . ${ssearch_path}
    """
}