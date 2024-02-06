# ONTViSc (ONT-based Viral Screening for Biosecurity)

## Introduction
eresearchqut/ontvisc is a Nextflow-based bioinformatics pipeline designed to help diagnostics of viruses and viroid pathogens for biosecurity. It takes fastq files generated from either amplicon or whole-genome sequencing using Oxford Nanopore Technologies as input.

The pipeline can either: 1) perform a direct search on the sequenced reads, 2) generate clusters, 3) assemble the reads to generate longer contigs or 4) directly map reads to a known reference. 

The reads can optionally be filtered from a plant host before performing downstream analysis.

**Sections:**

1. [Pipeline overview](#pipeline-overview)
2. [Running the pipeline](#running-the-pipeline)

## Pipeline overview
![diagram pipeline](docs/images/ONTViSc_pipeline.jpeg)

- Data quality check (QC) and preprocessing
  - Merge fastq files (optional)
  - Raw fastq file QC (Nanoplot)
  - Trim adaptors (PoreChop ABI - optional)
  - Filter reads based on length and/or quality (Chopper - optional)
  - Reformat fastq files so read names are trimmed after the first whitespace (bbmap)
  - Processed fastq file QC (if PoreChop and/or Chopper is run) (Nanoplot)
- Host read filtering
  - Align reads to host reference provided (Minimap2)
  - Extract reads that do not align for downstream analysis (seqtk)
- QC report
  - Derive read counts recovered pre and post data processing and post host filtering
- Read classification analysis mode
- Clustering mode
  - Read clustering (Rattle)
  - Convert fastq to fasta format (seqtk)
  - Cluster scaffolding (Cap3)
  - Megablast homology search against ncbi or custom database (blast)
  - Derive top candidate viral hits
- De novo assembly mode
  - De novo assembly (Canu or Flye)
  - Megablast homology search against ncbi or custom database or reference (blast)
  - Derive top candidate viral hits
- Read classification mode
  - Option 1 Nucleotide-based taxonomic classification of reads (Kraken2, Braken)
  - Option 2 Protein-based taxonomic classification of reads (Kaiju, Krona)
  - Option 3 Convert fastq to fasta format (seqtk) and perform direct homology search using megablast (blast)
- Map to reference mode
  - Align reads to reference fasta file (Minimap2) and derive bam file and alignment statistics (Samtools)

Detailed instructions can be found in [wiki](https://github.com/eresearchqut/ontvisc/wiki).
A step-by-step guide with instructions on how to set up and execute the ONTvisc pipeline on one of the HPC systems: Lyra (Queensland University of Technology), Setonix (Pawsey) and Gadi (National Computational Infrastructure) can be found [here](https://mantczakaus.github.io/ontvisc_guide/).


## Running the pipeline
- Run the command:
  ```
  nextflow run eresearchqut/ontvisc -profile {singularity, docker} --samplesheet index.csv
  ```
  The first time the command runs, it will download the pipeline into your assets.  

  The source code can also be downloaded directly from GitHub using the git command:
  ```
  git clone https://github.com/eresearchqut/ontvisc.git
  ```

- Provide an index.csv file.  
  Create a comma separated file that will be the input for the workflow. By default the pipeline will look for a file called “index.csv” in the base directory but you can specify any file name using the ```--samplesheet [filename]``` in the nextflow run command. This text file requires the following columns (which needs to be included as a header): ```sampleid,sample_files``` 

  **sampleid** will be the sample name that will be given to the files created by the pipeline  
  **sample_path** is the full path to the fastq files that the pipeline requires as starting input  

  This is an example of an index.csv file which specifies the name and path of fastq.gz files for 2 samples. Specify the full path length for samples with a single fastq.gz file. If there are multiple fastq.gz files per sample, place them all in a single folder and the path can be specified on one line using an asterisk:
  ```
  sampleid,sample_files
  MT212,/path_to_fastq_file_folder/*fastq.gz
  MT213,/path_to_fastq_file_folder/*fastq.gz
  ```

- Specify a profile:
  ```
  nextflow run eresearchqut/ontvisc -profile {singularity, docker} --samplesheet index_example.csv
  ```
  setting the profile parameter to one of ```docker``` or ```singularity``` to suit your environment.
  
- Specify one analysis mode: ```--analysis_mode {read classification, clustering, assembly, map2ref}``` (see below for more details)

- To set additional parameters, you can either include these in your nextflow run command:
  ```
  nextflow run eresearchqut/ontvisc -profile {singularity, docker} --samplesheet index_example.csv --adapter_trimming
  ```
  or set them to true in the nextflow.config file.
  ```
  params {
    adapter_trimming = true
  }
  ```

- A test is provided to check if the pipeline was successfully installed. The test.fastq.gz file is derived from of a plant infected with Miscanthus sinensis mosaic virus. To use the test, run the following command, selecting the adequate profile (singularity/docker):
  ```
  nextflow run eresearchqut/ontvisc -profile test,{singularity, docker}
  ```
The test requires 2 cpus at least 16Gb of memory to run and can be executed locally.  


The command should take one minute to run and nextflow should output the following log:
<p align="left"><img src="docs/images/nextflow_log.png"></p>

If the installation is successful, it will generate a results/test folder with the folloiwng structure:
```
results/
└── test
    ├── assembly
    │   ├── blast_to_ref
    │   │   └── blastn_reference_vs_flye_assembly.txt
    │   └── flye
    │       ├── test_flye_assembly.fasta
    │       ├── test_flye.fastq
    │       └── test_flye.log
    ├── preprocessing
    │   └── test_preprocessed.fastq.gz
    └── qc
        └── nanoplot
            └── test_raw_NanoPlot-report.html
``` 





### Authors
Marie-Emilie Gauthier <gauthiem@qut.edu.au>  
Craig Windell <c.windell@qut.edu.au>  
Magdalena Antczak <magdalena.antczak@qcif.edu.au>  
Roberto Barrero <roberto.barrero@qut.edu.au>  
