process GenerateSummaryReport {
    publishDir "${params.output_final_dir}", mode: "copy"


    input:
    path pseudogene_annotations_path
    val min_coverage

    output:
    path "summary_report.txt"

    script:
    """
    python3 ${workflow.projectDir}/../src/generate_summary_report.py ${pseudogene_annotations_path} ./summary_report.txt ${min_coverage}
    """
}