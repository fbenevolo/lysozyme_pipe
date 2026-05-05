process AnnotateRegionsWithBestProteins {
    publishDir "${params.output_final_dir}", mode: "copy"

    input:
    path genomic_regions_path
    path filtered_hits_path

    output:
    path "region_annotations.tsv"

    script:
    """
    python3 ${workflow.projectDir}/../src/score_density.py ${genomic_regions_path} ${filtered_hits_path} ./region_annotations.tsv 
    """
}