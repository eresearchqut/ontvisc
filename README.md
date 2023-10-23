# ONTViSc (ONT-based Viral Screening for Biosecurity)

## Introduction
eresearchqut/ontvisc is a Nextflow-based bioinformatics pipeline designed to help diagnostics of viruses and viroid pathogens for biosecurity. It takes fastq files generated from either amplicon or whole-genome sequencing using Oxford Nanopore Technologies as input.

The pipeline can either: 1) perform a direct search on the sequenced reads, 2) generate clusters, 3) assemble the reads to generate longer contigs or 4) directly map reads to a known reference. 

The reads can optionally be filtered from a plant host before performing downstream analysis.

## Pipeline overview
![diagram pipeline](docs/images/ONTViSc_pipeline.jpeg)

- Data quality check (QC) and preprocessing
  - NanoFilt - QC pre data processing
  - PoreChop ABI (optional)
  - Chopper (optional)
  - bbmap - Reformat fastq files so read names are trimmed after the first whitespace
  - NanoFilt - QC post data processing (if PoreChop and/or Chopper is run)
- Host read filtering
  - Minimap2 - align reads to host reference provided
  - seqtk - extract reads that do not align for downstream analysis
- QC report
  - custom python script derives read counts recovered pre and post data processing and post host filtering
- Read classification analysis mode
- Clustering mode
  - Rattle - clusters reads
  - seqtk - convert fastq to fasta format
  - Cap3 - scaffolds clusters
  - Blastn - perform megablast homology search against ncbi or custom database
  - custom python script - derive top candidate viral hits
- De novo assembly mode


### Authors
Marie-Emilie Gauthier <gauthiem@qut.edu.au>  
Craig Windell <c.windell@qut.edu.au>  
Magdalena Antczak <magdalena.antczak@qcif.edu.au>  
Roberto Barrero <roberto.barrero@qut.edu.au>  