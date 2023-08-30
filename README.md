# ONTViSc (ONT-based Viral Screening for Biosecurity)

## Introduction
eresearchqut/ontvisc is a Nextflow-based bioinformatics pipeline designed to help diagnostics of viruses and viroid pathogens for biosecurity. It takes fastq files generated from either amplicon or whole-genome sequencing using Oxford Nanopore Technologies as input.

The pipeline can either: 1) perform a direct search on the sequenced reads, 2) assemble the reads to generate longer contigs or 3) directly map reads to a known reference. 

The reads can optionally be filtered from a plant host before performing downstream analysis.

## Pipeline summary
![diagram pipeline](docs/images/OVISP_pipeline.jpeg)

## Run the Pipeline
1. Install Nextflow [`Nextflow`](https://www.nextflow.io/docs/latest/getstarted.html#installation)

2. Install [`Docker`](https://docs.docker.com/get-docker/) or [`Singularity`](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html#quick-installation-steps) to suit your environment.

3. Download the pipeline and test it on minimal datatests:
[ TO DO ]

4. Run with your own data

- Provide an index.csv file.  
  Create a TAB delimited text file that will be the input for the workflow. By default the pipeline will look for a file called “index.csv” in the base directory but you can specify any file name using the --indexfile [filename] in the nextflow run command. This text file requires the following columns (which needs to be included as a header): ```sampleid,sample_path``` 

  **sampleid** will be the sample name that will be given to the files created by the pipeline  
  **sample_path** is the full path to the fastq files that the pipeline requires as starting input  

  This is an example of an index.csv file which specifies the name and path of fastq.gz files for 2 samples. If there are multiple fastq.gz files in the folder, the path can be specified on one line using an asterisk:
  ```
  sampleid,sample_path
  MT212,/path_to_fastq_file/*fastq.gz
  MT213,/path_to_fastq_file/*fastq.gz
  ```

- Run the command:
  ```bash
  nextflow run main.nf -profile {singularity, docker} --indexfile index_example.csv
  ```
  setting the profile parameter to one of
  ```
  docker
  singularity
    ```  
  to suit your environment.

If you need to set additional parameters, you can either include these in your nextflow run command:
```
nextflow run main.nf -profile {singularity, docker} --indexfile index_example.csv --adapter_trimming
```

or set them to true in the nextflow.config file.
```
params {
  adapter_trimming = true
}
```

## Example of commands
1) running QC steps to have a preliminary look at the data before proceeding with downstream analysis.
```
nextflow run ~/path/to/ontvisc_repo/main.nf  --qc_only
```


2) remove adapters, perform de novo assembly with Canu and map the resulting contigs to a reference
```
nextflow run ~/path/to/ontvisc_repo/main.nf  -resume --adapter_trimming \\
                                                    --denovo_assembly --canu \\
                                                    --canu_options 'useGrid=false' \\
                                                    --canu_genome_size [genome size of virus target] \\
                                                    --blast_vs_ref  \\
                                                    --reference /path/to/reference/reference.fasta
```