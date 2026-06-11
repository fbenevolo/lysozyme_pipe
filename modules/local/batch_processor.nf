process BatchProcessor {
    publishDir "${params.output_dir}", mode: "copy"
    
    input:
    path input_dir
    path lysozymes_fasta
    val output_dir
    val min_identity
    val min_score
    val min_disablements
    val min_coverage
    val final_min_identity
    val num_threads

    output:
    path "${params.output_dir}"

    script:
    """
    python3 ${workflow.projectDir}/../src/batch_processor.py \
        ${input_dir} \
        ${lysozymes_fasta} \
        ${output_dir} \
        ${min_identity} \
        ${min_score} \
        ${min_disablements} \
        ${min_coverage} \
        ${final_min_identity} \
        ${num_threads}
    """
}