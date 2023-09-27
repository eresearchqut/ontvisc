# ONTViSc (ONT-based Viral Screening for Biosecurity)

###Authors
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

| Mode | Index | Description |
| --- | --- | --- |
| host_filtering | host_fasta | reference fasta file |
| clustering blast_mode localdb | blastn_db | viral database e.g. [`RVDB`](https://rvdb.dbi.udel.edu/) |
| denovo_assembly blast_mode localdb | | |
| read_classification  blast_mode localdb | | |
| clustering blast_mode ncbi | blastn_db | NCBI NT |
| denovo_assembly blast_mode ncbi | | |
| read_classification  blast_mode ncbi | | |
| read_classification kraken2 | krkdb | kraken index |
| read_classification kaiju | kaiju_dbname | path to kaiju_db_*.fmi |
|                           | kaiju_nodes | path to nodes.dmp |
|                           | kaiju_names | path to names.dmp |
| map2ref or blast_vs_ref | reference | viral reference sequence |

- If you have access to a host genome reference or sequences and want to filter your reads against it/them before running your analysis, you will have to specify the ``--host_filtering``parameter and provide the path to the host fasta file with ``--host_fasta /path/to/host/fasta/file``

- If you  want to run homology searches against a viral database (e.g. [`RVDB`](https://rvdb.dbi.udel.edu/), you will need to specify the ``--blast_mode localdb`` parameter and provide the path to the database by specifying: ``--blastn_db /path/to/viral/db``. 

- If you  want to run homology searches against the public NCBI NT database instead, you need to set the parameter ```--blast_mode ncbi```
This parameter is set by default in the nextflow.config file:
```
params {
  blast_mode = ncbi
}
```

Download the NCBI NT database locally, following the detailed steps available at https://www.ncbi.nlm.nih.gov/books/NBK569850/. Create a folder where you will store your NCBI databases. It is good practice to include the date of download. For instance:
```
mkdir blastDB/20230930
```
You will need to use the update_blastdb.pl script from the blast+ version used with the pipeline.
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

- Run the command:
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
- Run only the quality control step to have a preliminary look at the data before proceeding with downstream analysis.
```
nextflow run ~/path/to/ontvisc_repo/main.nf  --qc_only
```

- Trim adapters using [`PoreChop ABI`](https://github.com/rrwick/Porechop)
```
nextflow run ~/path/to/ontvisc_repo/main.nf  --adapter_trimming
```

Additional PoreChop parameters can be specified using ```--porechop_options '{options}'```. Please refer to PoreChop manual.

- If the data analysed was derived using RACE reactions, a final primer check can be performed after the de novo assembly step using the ```--final_primer_check``` option. The pipeline will check for the presence of any residual universal RACE primers at the end of the assembled contigs.

- Perform a quality filtering step  using ```--qual_filt``` and either the ```chopper``` (default) or the ```nanoFilt``` option.
Additional Chopper and NanoFilt parameters can be specified using ```--chopper_options``` and ```--nanofilt_options``` respectively.

```
nextflow run ~/path/to/ontvisc_repo/main.nf  --qual_filt
```

5. Provide required databases




## Example of commands


1) check for presence of adapters, perform de novo assembly with Canu and map the resulting contigs to a reference.
If you do not know the size of your targetted genome, you can ommit the ```--canu_genome_size parameter```. However, if your sample is likely to contain a lot of plant RNA/DNA material, we recommend providing an approximate genome size. For instance RNA viruses are on average 10 kb in size (see [`Holmes 2009`](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC2954018/))

```
nextflow run ~/path/to/ontvisc_repo/main.nf  -resume --adapter_trimming \
                                                    --denovo_assembly --canu \
                                                    --canu_options 'useGrid=false' \
                                                    --canu_genome_size [genome size of virus target] \
                                                    --blast_vs_ref  \
                                                    --reference /path/to/reference/reference.fasta
```

2) check for presence of adapters and perform a direct read homology search using megablast and the NCBI NT database. You will need to download a local copy of the NCBI NT database. The blast search will be split into several jobs, containing 10,000 reads each, that will run in parallel. The pipeline will use 8 cpus when running the blast process.

```
nextflow run ~/path/to/ontvisc_repo/main.nf  -resume --adapter_trimming \
                                                    --read_classification \
                                                    --megablast \
                                                    --blast_threads 8 \
                                                    --blast_mode ncbi \
                                                    --blastn_db /path/to/ncbi_blast_db/nt
```

3) check for presence of adapters and this time perform a direct read homology search using megablast against a viral database.

```
nextflow run ~/path/to/ontvisc_repo/main.nf  -resume --adapter_trimming \
                                                    --read_classification \
                                                    --megablast \
                                                    --blast_threads 8 \
                                                    --blast_mode localdb \
                                                    --blastn_db /path/to/local_blast_db
```

4) check for presence of adapters and perform a direct taxonomic classification of reads using Kraken2 and Kaiju. You will need to download Kraken2 index (e.g. PlusPFP) and Kaiju indexes (e.g. kaiju_db_rvdb).


```
nextflow run ~/path/to/ontvisc_repo/main.nf  -resume --adapter_trimming \
                                                    --read_classification \
                                                    --kraken2 \
                                                    --krkdb = /path/to/kraken2_db \
                                                    --kaiju \
                                                    --kaiju_dbname /path/to/kaiju/kaiju.fmi \
                                                    --kaiju_nodes /path/to/kaiju/nodes.dmp \
                                                    --kaiju_names /path/to/kaiju/names.dmp
```
