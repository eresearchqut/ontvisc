# OViSP (ONT-based Viral Screening for Plant diagnostics)

## Introduction
eresearchqut/OViSP is a Nextflow-based bioinformatics pipeline. It was designed to help phytosanitary diagnostics of viruses and viroid pathogens. It takes fastq files generated from amplicon or whole-genome sequencing using Oxford Nanopore Technologies as input.

The pipeline can either perform a direct homology search on the sequenced reads or assemble the reads to generate longer contigs before annotation and diagnosis. The reads can optionally be filtered from a plant host before performing the homology search.

## Pipeline summary
![diagram pipeline](docs/images/OVISP_pipeline.jpg)