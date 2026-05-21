process SavePseudogeneAnnotations {
    publishDir "${params.output_final_dir}", mode: "copy"

    input:
    path pseusogene_annotations_path

    output:
    path "pseudogene_annotations_final.tsv"
    
    script:
    """
    python3 ${workflow.projectDir}/../src/save_pseudogene_annotations.py ${pseusogene_annotations_path} ./pseudogene_annotations_final.tsv
    """
}