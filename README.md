## EthSEQ Reference Model Generation Pipeline

Snakemake pipeline for the creation of [EthSEQ](https://github.com/cibiobcg/EthSEQ) reference models. Automate model generation for various subpopulations and WES kits, with parallelization for efficiency.

- Flexible configuration (path.yaml)
- Model creation for diverse datasets
- Parallelization for faster processing
- Not mandatory for EthSEQ model creation (alternative methods available into EthSEQ)

## Required Input Files

To use this pipeline, you'll need the following files:

- BED Files: WES kit region files in BED format. These define the genomic regions targeted by your Whole Exome Sequencing (WES) kit.
- VCF Files: Variant Call Format (VCF) files containing genotype data for individuals with known ethnicities. These are your training data for the model.
- Sample Information File: A file (e.g., CSV or TSV) containing sample identifiers and their corresponding ethnicities. This file links the genotype data in the VCF files to the known ethnicities.
- dbSNPs dataset file: A file containing variant information, to add rsID to the model
