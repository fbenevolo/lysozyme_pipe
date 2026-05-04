process BedtoolsSaveMergedRegions {
    publishDir "${params.output_merge}", mode: "copy"

    input:
    path genomic_regions_path
    val merged_regions_path
    val genome_id

    output:
    path "merged_regions.tsv"

    script:
    """
    python3 ${workflow.projectDir}/../src/bedtools_save_merged_regions.py ${genomic_regions_path} ${merged_regions_path} ${genome_id}
    """
}