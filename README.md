# Nanopore Sequencing Analysis Pipeline

This repository contains a Nextflow-based bioinformatics pipeline for the analysis of Nanopore sequencing data.

## Features

- Adapter trimming using Porechop
- Alignment to a reference genome with Minimap2
- Variant calling with Bcftools
- Consensus sequence generation with iVar

## Dependencies

The pipeline is dependent on the following software:

- Nextflow
- Docker

The required bioinformatics tools (Porechop, Minimap2, Samtools, Bcftools, iVar) are encapsulated within Docker containers, which are automatically fetched at runtime.

## Usage

Ensure that both Nextflow and Docker are installed on your system. 

To run the pipeline, use the following command:

```bash
nextflow run main.nf --reads '/path/to/reads/*.fastq' --ref '/path/to/reference.fasta' --outdir 'output_directory' -profile docker
