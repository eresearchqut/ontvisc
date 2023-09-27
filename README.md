# ONTViSc (ONT-based Viral Screening for Biosecurity)

### Authors
Marie-Emilie Gauthier <gauthiem@qut.edu.au>

## Introduction
eresearchqut/ontvisc is a Nextflow-based bioinformatics pipeline designed to help diagnostics of viruses and viroid pathogens for biosecurity. It takes fastq files generated from either amplicon or whole-genome sequencing using Oxford Nanopore Technologies as input.

The pipeline can either: 1) perform a direct search on the sequenced reads, 2) assemble the reads to generate longer contigs or 3) directly map reads to a known reference. 

The reads can optionally be filtered from a plant host before performing downstream analysis.

## Pipeline summary
![diagram pipeline](docs/images/ONTViSc_pipeline.jpeg)

## Installation
1. Install Nextflow [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. Install [`Docker`](https://docs.docker.com/get-docker/) or [`Singularity`](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps) to suit your environment.

3. Download the pipeline and test it on minimal datatests:
The source code can be downloaded directly from GitHub using the git command line client:
```
git clone https://github.com/maelyg/ontvisc.git
```

## Installing the required indexes and references
Depending on the mode you are intersted to run, you will need to install some databases and references.

| Mode | Index | Description |
| --- | --- | --- |
| --host_filtering | --host_fasta | path to host fasta file to use for read filtering|
| --blast_vs_ref | --reference | path to viral reference sequence fasta file to perform homology search on reads (read_classification), clusters (clustering) or contigs (de novo) |
| --blast_mode localdb | --blastn_db | path to viral blast database e.g. [`RVDB`](https://rvdb.dbi.udel.edu/) to perform homology search on reads (read_classification), clusters (clustering) or contigs (de novo)|
| --blast_mode ncbi | --blastn_db | path to NCBI nt database, taxdb.btd and taxdb.bti to perform homology search on reads (read_classification), clusters (clustering) or contigs (de novo)|
| --read_classification --kraken2 | --krkdb | path to kraken index folder e.g. PlusPFP|
| --read_classification --kaiju | --kaiju_dbname | path to kaiju_db_*.fmi |
|                           | --kaiju_nodes | path to nodes.dmp |
|                           | --kaiju_names | path to names.dmp |
| --map2ref | --reference | path to viral reference sequence fasta file to perform alignment |

- If you have access to a host genome reference or sequences and want to filter your reads against it/them before running your analysis, specify the ``--host_filtering`` parameter and provide the path to the host fasta file with ``--host_fasta /path/to/host/fasta/file``.

- If you want to run homology searches against a viral database (e.g. [`RVDB`](https://rvdb.dbi.udel.edu/), you will need to specify the ``--blast_mode localdb`` parameter and provide the path to the database by specifying: ``--blastn_db /path/to/viral/db``. 

- If you want to run homology searches against the public NCBI NT database instead, you need to set the parameter ```--blast_mode ncbi```
This parameter is set by default in the nextflow.config file:
```
params {
  blast_mode = ncbi
}
```

Download a local copy of the NCBI NT database, following the detailed steps available at https://www.ncbi.nlm.nih.gov/books/NBK569850/. Create a folder where you will store your NCBI databases. It is good practice to include the date of download. For instance:
```
mkdir blastDB/20230930
```
You will need to use a current update_blastdb.pl script from the blast+  version used with the pipeline (ie 2.13.0).
For example:
```
perl update_blastdb.pl --decompress nt [*]
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

- Provide an index.csv file.  
  Create a TAB delimited text file that will be the input for the workflow. By default the pipeline will look for a file called “index.csv” in the base directory but you can specify any file name using the ```--samplesheet [filename]``` in the nextflow run command. This text file requires the following columns (which needs to be included as a header): ```sampleid,sample_files``` 

  **sampleid** will be the sample name that will be given to the files created by the pipeline  
  **sample_path** is the full path to the fastq files that the pipeline requires as starting input  

  This is an example of an index.csv file which specifies the name and path of fastq.gz files for 2 samples. If there are multiple fastq.gz files in the folder, the path can be specified on one line using an asterisk:
  ```
  sampleid,sample_files
  MT212,/path_to_fastq_file/*fastq.gz
  MT213,/path_to_fastq_file/*fastq.gz
  ```

- Specify a profile:
  ```bash
  nextflow run main.nf -profile {singularity, docker} --samplesheet index_example.csv
  ```
  setting the profile parameter to one of
  ```
  docker
  singularity
    ```  
  to suit your environment.

To set additional parameters, you can either include these in your nextflow run command:
```
nextflow run main.nf -profile {singularity, docker} --samplesheet index_example.csv --adapter_trimming
```

or set them to true in the nextflow.config file.
```
params {
  adapter_trimming = true
}
```

# Running QC step
By default the pipeline will run a quality control check of the raw reads using NanoPlot.

- Run only the quality control step to have a preliminary look at the data before proceeding with downstream analysis by specifying the ```--qc_only``` parameter.

# Preparing reads before analysis
Raw read can be trimmed of adapters and/or quality filtered.
- Trim adapters can be trimmed using [`PoreChop ABI`](https://github.com/rrwick/Porechop) by specifying the ``--adapter_trimming`` parameter. PoreChop ABI options can be specified using ```--porechop_options '{options}'```. Please refer to the PoreChop manual.

- If the data analysed was derived using RACE reactions, a final primer check can be performed after the de novo assembly step using the ```--final_primer_check``` option. The pipeline will check for the presence of any residual universal RACE primers at the end of the assembled contigs.

- Perform a quality filtering step using ```--qual_filt``` and run either [`Chopper`](https://github.com/wdecoster/chopper) or [`NanoFilt`](https://github.com/wdecoster/nanofilt) by specifying the ```chopper``` (default) or the ```nanofilt``` option respectively. Chopper and NanoFilt options can be specified using the ```--chopper_options``` and the ```--nanofilt_options``` respectively. Please refer to the Chopper and NanoFilt manuals.

- If you trim raw read of adapters and/or quality filter the raw reads, an additional quality control step will be performed and a qc report will be generated summarising the read counts retained at each step.

# Filtering host reads
- Reads mapping to a host genome reference or sequences can be filtered out by specifying the ``--host_filtering`` parameter and provide the path to the host fasta file with ``--host_fasta /path/to/host/fasta/file``.

# Running read classification (--read_classification)
- Perform a direct blast homology search using megablast (``--megablast``).

Example 1 using a viral database:
```
# Check for presence of adapters.
# Perform a direct read homology search using megablast against a viral database.

nextflow run ~/path/to/ontvisc_repo/main.nf -resume --adapter_trimming \
                                                    --read_classification \
                                                    --megablast \
                                                    --blast_threads 8 \
                                                    --blast_mode localdb \
                                                    --blastn_db /path/to/local_blast_db
```

Example 2 using NCBI nt:
```
# Check for presence of adapters .
# Perform a direct read homology search using megablast and the NCBI NT database. 
# You will need to download a local copy of the NCBI NT database. 
# The blast search will be split into several jobs, containing 10,000 reads each, that will run in parallel. 
# The pipeline will use 8 cpus when running the blast process.

nextflow run ~/path/to/ontvisc_repo/main.nf -resume --adapter_trimming \
                                                    --read_classification \
                                                    --megablast \
                                                    --blast_threads 8 \
                                                    --blast_mode ncbi \
                                                    --blastn_db /path/to/ncbi_blast_db/nt
```

- Perform a direct taxonomic classification of reads using Kraken2 and Kaiju
Example:
```
# Check for presence of adapters
# Perform a direct taxonomic classification of reads using Kraken2 and Kaiju. 
# You will need to download Kraken2 index (e.g. PlusPFP) and Kaiju indexes (e.g. kaiju_db_rvdb).

nextflow run ~/path/to/ontvisc_repo/main.nf -resume --adapter_trimming \
                                                    --read_classification \
                                                    --kraken2 \
                                                    --krkdb = /path/to/kraken2_db \
                                                    --kaiju \
                                                    --kaiju_dbname /path/to/kaiju/kaiju.fmi \
                                                    --kaiju_nodes /path/to/kaiju/nodes.dmp \
                                                    --kaiju_names /path/to/kaiju/names.dmp
```

# Run de novo assembly (--denovo_assembly)
You can run a de novo assembly using either either Flye or Canu. 

- Canu:
Canu options can be specified using the ```--canu_options``` parameter.
If you do not know the size of your targetted genome, you can ommit the ```--canu_genome_size [genome size of virus target]```. However, if your sample is likely to contain a lot of plant RNA/DNA material, we recommend providing an approximate genome size. For instance RNA viruses are on average 10 kb in size (see [`Holmes 2009`](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2954018/))

You can perform a homology search against the contigs generated using either a viral genome reference, a viral database or NCBI nt.

Example:
```
# Check for presence of adapters
# Perform de novo assembly with Canu
# Blast the resulting contigs to a reference.
nextflow run ~/path/to/ontvisc_repo/main.nf -resume --adapter_trimming \
                                                    --denovo_assembly --canu \
                                                    --canu_options 'useGrid=false' \
                                                    --canu_genome_size 0.01m \
                                                    --blast_vs_ref  \
                                                    --reference /path/to/reference/reference.fasta
```

# Run clustering (--clustering)
Run the clustering tool [`RATTLE`](https://github.com/comprna/RATTLE#Description-of-clustering-parameters)
Set the parameter ```--rattle_clustering_options '--raw'``` to use all the reads without any length filtering during the RATTLE clustering step if your amplicon is known to be shorter than 150 bp.

When the amplicon is of known size, we recommend setting up the parameters ```--lower-length [number]``` (default: 150) and ```--upper-length [number]``` (default: 100,000) to filter out reads shorter and longer than the expected size.

Set the parameter ```--rattle_clustering_options '--rna'``` and ```--rattle_polishing_options '--rna'``` if the data is direct RNA (disables checking both strands).

Example:
```
# With this command all reads will be retained during the clustering step.
nextflow run ~/path/to/ontvisc_repo/main.nf -resume --clustering \
                                                     --rattle_clustering_options '--raw' \
                                                     --blast_threads 8 \
                                                     --blastn_db /path/to/ncbi_blast_db/nt
```