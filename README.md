# ONTViSc (ONT-based Viral Screening for Biosecurity)

## Introduction
eresearchqut/ontvisc is a Nextflow-based bioinformatics pipeline designed to help diagnostics of viruses and viroid pathogens for biosecurity. It takes fastq files generated from either amplicon or whole-genome sequencing using Oxford Nanopore Technologies as input.

The pipeline can either: 1) perform a direct search on the sequenced reads, 2) generate clusters, 3) assemble the reads to generate longer contigs or 4) directly map reads to a known reference. 

The reads can optionally be filtered from a plant host before performing downstream analysis.

**Sections:**

1. [Pipeline overview](#pipeline-overview)  
2. [Installation](#installation)  
a. [Requirements](#requirements)  
b. [Installing the required indexes and references](#installing-the-required-indexes-and-references)  
3. [Running the pipeline](#running-the-pipeline)  
a. [Run test data](#run-test-data)  
b. [QC step](#qc-step)  
c. [Preprocessing reads](#preprocessing-reads)  
d. [Host read filtering](#host-read-filtering)  
5. [Authors](#authors)

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

## Installation
### Requirements  
1. Install Nextflow [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. Install [`Docker`](https://docs.docker.com/get-docker/) or [`Singularity`](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps) to suit your environment.

### Installing the required indexes and references
Depending on the mode you are interested to run, you will need to install some databases and references.

| Mode | Index | Description |
| --- | --- | --- |
| --host_filtering | --host_fasta | path to host fasta file to use for read filtering|
| --blast_vs_ref | --reference | path to viral reference sequence fasta file to perform homology search on reads (read_classification), clusters (clustering) or contigs (de novo) |
| --blast_mode localdb | --blastn_db | path to [`viral blast database`](https://zenodo.org/records/10117282) to perform homology search on reads (read_classification), clusters (clustering) or contigs (de novo)|
| --blast_mode ncbi | --blastn_db | path to NCBI nt database, taxdb.btd and taxdb.bti to perform homology search on reads (read_classification), clusters (clustering) or contigs (de novo)|
| --read_classification --kraken2 | --krkdb | path to kraken index folder e.g. PlusPFP|
| --read_classification --kaiju | --kaiju_dbname | path to kaiju_db_*.fmi |
|                           | --kaiju_nodes | path to nodes.dmp |
|                           | --kaiju_names | path to names.dmp |
| --map2ref | --reference | path to viral reference sequence fasta file to perform alignment |

- If you have access to a host genome reference or sequences and want to filter your reads against it/them before running your analysis, specify the `--host_filtering` parameter and provide the path to the host fasta file with `--host_fasta /path/to/host/fasta/file`.

- The homology searches is set by default against the public NCBI NT database in the nextflow.config file (`--blast_mode ncbi`)

<<<<<<< HEAD
- If you want to run homology searches against a viral database instead, you will need to download it [`here`](https://zenodo.org/records/10183620) by using the following steps:  
```
wget https://zenodo.org/records/10183620/files/VirDB_20230913.tar.gz?download=1
tar -xf VirDB_20230913.fasta.tar.gz
```
Specify the ``--blast_mode localdb`` parameter and provide the path to the database by specifying ``--blastn_db /path/to/viral/db``.

  Download a local copy of the NCBI NT database, following the detailed steps available at https://www.ncbi.nlm.nih.gov/books/NBK569850/. Create a folder where you will store your NCBI databases. It is good practice to include the date of download. For instance:
  ```
  mkdir blastDB/20230930
  ```
  You will need to use a current update_blastdb.pl script from the blast+  version used with the pipeline (ie 2.13.0).
  For example:
  ```
  perl update_blastdb.pl --decompress nt
  perl update_blastdb.pl taxdb
  tar -xzf taxdb.tar.gz
  ```

  Make sure the taxdb.btd and the taxdb.bti files are present in the same directory as your blast databases.
  Specify the path of your local NCBI blast nt directories in the nextflow.config file.
  For instance:
  ```
  params {
    --blastn_db = '/work/hia_mt18005_db/blastDB/20230930/nt'
  }
  ```
- To run nucleotide taxonomic classification of reads using Kraken2, download the pre-built index relevant to your data and provided by [`Kraken2`](https://benlangmead.github.io/aws-indexes/k2) (for example, PlusPFP can be chosen for searching viruses in plant samples).  

- To run protein taxonomic classification using Kaiju, download the pre-built index relevant to your data. Indexes are listed on the README page of [`Kaiju`](https://github.com/bioinformatics-centre/kaiju) (for example refseq, refseq_nr, refseq_ref, progenomes, viruses, nr, nr_euk or rvdb). After the download is finished, you should have 3 files: kaiju_db_*.fmi, nodes.dmp, and names.dmp, which are all needed to run Kaiju.
You will have to specify the path to each of these files (using the ``--kaiju_dbname``, the ``--kaiju_nodes`` and the ``--kaiju_names`` parameters respectively.

- If you want to align your reads to a reference genome (--map2ref) or blast against a reference (--blast_vs_ref), you will have to specify its path using `--reference`.


## Running the pipeline
### Run test data
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

If the installation is successful, it will generate a results/test folder with the following structure:
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

## QC step
By default the pipeline will run a quality control check of the raw reads using NanoPlot.

- Run only the quality control step to have a preliminary look at the data before proceeding with downstream analyses by specifying the ```--qc_only``` parameter.

## Preprocessing reads
If multiple fastq files exist for a single sample, they will first need to be merged using the `--merge` option.
Then the read names of the fastq file created will be trimmed after the first whitespace, for compatiblity purposes with all downstream tools.  

Reads can also be optionally trimmed of adapters and/or quality filtered:  
- Search for presence of adapters in sequences reads using [`Porechop ABI`](https://github.com/rrwick/Porechop) by specifying the ``--adapter_trimming`` parameter. Porechop ABI parameters can be specified using ```--porechop_options '{options} '```, making sure you leave a space at the end before the closing quote. Please refer to the Porechop manual.  
To limit the search to known adapters listed in [`adapter.py`](https://github.com/bonsai-team/Porechop_ABI/blob/master/porechop_abi/adapters.py), just specify the ```--adapter_trimming``` option.  
To search ab initio for adapters on top of known adapters, specify ```--adapter_trimming --porechop_options '-abi '```.  
To limit the search to custom adapters, specify ```--adapter_trimming --porechop_custom_primers --porechop_options '-ddb '``` and list the custom adapters in the text file located under bin/adapters.txt following the format:
    ```
     line 1: Adapter name
     line 2: Start adapter sequence
     line 3: End adapter sequence
     --- repeat for each adapter pair---
     ```

- Perform a quality filtering step using [`Chopper`](https://github.com/wdecoster/chopper) by specifying the ```--qual_filt``` parameter. Chopper parameters can be specified using the ```--chopper_options '{options}'```. Please refer to the Chopper manual.  
For instance to filter reads shorter than 1000 bp and longer than 20000 bp, and reads with a minimum Phred average quality score of 10, you would specify: ```--qual_filt --chopper_options '-q 10 -l 1000 --maxlength 20000'```.  

A zipped copy of the resulting preprocessed and/or quality filtered fastq file will be saved in the preprocessing folder.  

If you trim raw read of adapters and/or quality filter the raw reads, an additional quality control step will be performed and a qc report will be generated summarising the read counts recovered before and after preprocessing for all samples listed in the index.csv file.

## Host read filtering
- Reads mapping to a host genome reference or sequences can be filtered out by specifying the ``--host_filtering`` parameter and provide the path to the host fasta file with ``--host_fasta /path/to/host/fasta/file``.

A qc report will be generated summarising the read counts recovered after host filtering.

## Authors
Marie-Emilie Gauthier <gauthiem@qut.edu.au>  
Craig Windell <c.windell@qut.edu.au>  
Magdalena Antczak <magdalena.antczak@qcif.edu.au>  
Roberto Barrero <roberto.barrero@qut.edu.au>  
