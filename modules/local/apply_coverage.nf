process ApplyCoverage {
    publishDir "${params.output_final_dir}", mode: "copy"

    input:
    path pseudogene_annotations_path
    val min_coverage

    output:
    path "pseudogene_annotations_with_coverage.jsonl"
    
    script:
    """
    python3 ${workflow.projectDir}/../src/apply_coverage_filter.py ${pseudogene_annotations_path} ${min_coverage} ./pseudogene_annotations_with_coverage.jsonl
    """
}