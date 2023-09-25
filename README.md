# ONTViSc (ONT-based Viral Screening for Biosecurity)

## Introduction
eresearchqut/ontvisc is a Nextflow-based bioinformatics pipeline designed to help diagnostics of viruses and viroid pathogens for biosecurity. It takes fastq files generated from either amplicon or whole-genome sequencing using Oxford Nanopore Technologies as input.

The pipeline can either: 1) perform a direct search on the sequenced reads, 2) assemble the reads to generate longer contigs or 3) directly map reads to a known reference. 

The reads can optionally be filtered from a plant host before performing downstream analysis.

## Pipeline summary
![diagram pipeline](docs/images/ONTViSc_pipeline.jpeg)

## Run the Pipeline
1. Install Nextflow [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. Install [`Docker`](https://docs.docker.com/get-docker/) or [`Singularity`](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps) to suit your environment.

3. Download the pipeline and test it on minimal datatests:
[ TO DO ]

4. Run with your own data

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

If you need to set additional parameters, you can either include these in your nextflow run command:
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


- If the data analysed used RACE reactions, a final primer check can be performed after de novo assembly using the ```--final_primer_check``` option. The pipeline will check for the presence of any residual universal RACE primers at the end of the assembled contigs.

- Provide a database
If you also want to run homology searches against public NCBI databases, you need to set the parameter ```--blast_mode ncbi```
This parameter is set by default in the nextflow.config file:
```
params {
  blast_mode = ncbi
}
```

Download these locally, following the detailed steps available at https://www.ncbi.nlm.nih.gov/books/NBK569850/. Create a folder where you will store your NCBI databases. It is good practice to include the date of download. For instance:
```
mkdir blastDB/20231130
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
  blast_db_dir = '/work/hia_mt18005_db/blastDB/20231130'
}
```

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

2) check for presence of adapters and perform a direct read homology search using megablast. You will need to download a local copy of the NCBI NT database. The blast search will be split into several jobs, containing 10,000 reads each, that will run in parallel. The pipeline will use 8 cpus when running the blast process.

```
nextflow run ~/path/to/ontvisc_repo/main.nf  -resume --adapter_trimming \
                                                    --read_classification \
                                                    --megablast \
                                                    --blast_threads 8 --blast_mode ncbi \
                                                    --blastn_db /path/to/ncbi_blast_db/nt \
```
