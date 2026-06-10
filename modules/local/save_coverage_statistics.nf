process SaveCoverageStatistics {
    publishDir "${params.output_final_dir}", mode: "copy"

    input:
    path pseudogene_annotations_path

    output: 
    path "coverage_statistics.csv"

    script:
    """
    python3 ${workflow.projectDir}/../src/save_coverage_statistics.py ${pseudogene_annotations_path} ./coverage_statistics.csv
    """
}