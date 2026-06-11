# Lysozyme Annotation Pipeline - Example Data

This directory contains minimal example data for testing the pipeline.

## Files Included

- `example_genome.fasta` - Small bacterial genome segment (~100kb)
- `reference_lysozymes.fasta` - Curated lysozyme protein sequences
- `example_metadata.tsv` - Metadata template for comparative analysis

## Quick Start

```bash
# Single genome analysis
python pipeline.py \
    -g examples/example_genome.fasta \
    -l examples/reference_lysozymes.fasta \
    -o examples/output/ \
    --num-threads 4

# View results
cat examples/output/final/pseudogene_annotations.tsv
```

## Getting Real Data

For real analysis, download genomes from:
- NCBI GenBank: https://www.ncbi.nlm.nih.gov/genome/
- ENA: https://www.ebi.ac.uk/ena/

Download lysozyme sequences from:
- UniProt: https://www.uniprot.org/ (search "lysozyme reviewed:yes")
