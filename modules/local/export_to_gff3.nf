process ExportToGFF3 {
    publishDir "${params.output_final_dir}", mode: "copy"

    input:
    path pseudogene_annotations_path
    val genome_id

    script:
    """
    python3 ${workflow.projectDir}/../src/export_gff3.py ${pseudogene_annotations_path} ${genome_id} ./lysozyme_annotations.gff3 
    """
}