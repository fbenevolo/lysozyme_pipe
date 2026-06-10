process AnnotatePseudogenes {
    publishDir "${params.output_final_dir}", mode: "copy"

    input:
    path regions_with_best_proteins
    path genome_fasta
    val min_disablements

    output:
    path "initial_pseudogene_annotation.jsonl"
    
    script:
    """
    python3 ${workflow.projectDir}/../src/pseudogene_detection.py ${regions_with_best_proteins} ${genome_fasta} ${min_disablements} "initial_pseudogene_annotation.jsonl" 
    """
}