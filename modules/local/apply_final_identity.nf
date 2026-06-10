process ApplyFinalIdentity {
    publishDir "${params.output_final_dir}", mode: "copy"
    
    input:
    path pseudogene_annotation_path
    val final_min_identity

    output: 
    path "pseudogene_annotations_with_final_identity.jsonl"

    script:
    """
    python3 ${workflow.projectDir}/../src/apply_final_identity_filter.py ${pseudogene_annotation_path} ${final_min_identity} ./pseudogene_annotations_with_final_identity.jsonl
    """
}