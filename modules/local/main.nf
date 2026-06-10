process Main {
    debug true

    script:
    """
    python3 ${projectDir}/../pipeline.py --genome ${projectDir}/../${params.genome_fasta} \
    --lysozymes ${projectDir}/../${params.lysozymes_fasta} --output ${params.output_folder}
    """
}