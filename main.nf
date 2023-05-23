#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage () {
    log.info """
    OViSP - ONT-based Viral Screening for Plants
    Marie-Emilie Gauthier 23/05/2023

    Usage:
    Run the command
    nextflow run eresearchqut/ovisp {optional arguments}...

    Optional arguments:
      -resume                           Resume a failed run
      --outdir                          Path to save the output file
                                        'results'
      --samplesheet '[path/to/file]'    Path to the csv file that contains the list of
                                        samples to be analysed by this pipeline.
                              Default:  'index.csv'
      Contents of samplesheet csv:
        sampleid,sample_files,reference
        SAMPLE01,/user/folder/*.fastq.gz,/path/to/reference.fasta

        sample_files can refer to a folder with a number of
        files that will be merged in the pipeline

        --flye_read_error               adjust parameters for given read error rate (as fraction e.g. 0.03)
                              Default:  0.03

        --flye_ont_mode                 Select from nano-raw, nano-corr, nano-hq
                              Default:  'nano-hq'

        --nanoq_code_start              Start codon position in the reference sequence
                              Default:  1

        --nanoq_read_length             Length cut off for read size
                              Default:  9000

        --nanoq_num_ref                 Number of references used in the alignment
                              Default:  1

        --nanoq_qual_threshhold         Base quality score cut off
                              Default:  5

        --nanoq_jump                    Increase this to make larger read intervals
                              Default:  10

    """.stripIndent()
}
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}


process MERGE {
  //publishDir "${params.outdir}/${sampleid}/merge", pattern: '*.fastq.gz', mode: 'link'
  tag "${sampleid}"
  label 'small'

  input:
    tuple val(sampleid), path(lanes), path(reference)
  output:
    tuple val(sampleid), path("${sampleid}.fastq.gz"), path(reference), emit: merged
  script:
  """
  cat ${lanes} > ${sampleid}.fastq.gz
  """

}

process NANOPLOT {
  publishDir "${params.outdir}/${sampleid}/nanoplot",  pattern: '*.html', mode: 'link', saveAs: { filename -> "${sampleid}_$filename" }
  publishDir "${params.outdir}/${sampleid}/nanoplot",  pattern: '*.NanoStats.txt', mode: 'link', saveAs: { filename -> "${sampleid}_$filename" }
  tag "${sampleid}"
  cpus 2

  container 'quay.io/biocontainers/nanoplot:1.41.0--pyhdfd78af_0'
  
  input:
    tuple val(sampleid), path(sample), path(reference)
  output:
    path("*.html")
    path("*NanoStats.txt")
  script:
  """
  NanoPlot -t 2 --fastq ${sample} --prefix raw --plots dot --N50
  """
}

process CUTADAPT_RACE {
  //publishDir "${params.outdir}/${sampleid}/canu", pattern: '*_cutadapt_filtered.fastq.gz', mode: 'link'
  publishDir "${params.outdir}/${sampleid}/canu", pattern: '*_cutadapt.log', mode: 'link'
  tag "${sampleid}"
  label 'medium'

  container 'quay.io/biocontainers/cutadapt:4.1--py310h1425a21_1'

  input:
    tuple val(sampleid), path(sample), path(reference)
  output:
    //path("${sampleid}_cutadapt_filtered.fastq.gz")
    path("${sampleid}_cutadapt.log")
    tuple val(sampleid), path("${sampleid}_cutadapt_filtered.fastq.gz"), path(reference), emit: cutadapt_filtered
  script:
  """
  if [[ ${params.race5} == true ]]; then
    cutadapt -j ${task.cpus} --times 3 -e 0.2 -g "AAGCAGTGGTATCAACGCAGAGTACGCGGG;min_overlap=14" -o ${sampleid}_filtered1.fastq.gz ${sample} > ${sampleid}_cutadapt.log
    cutadapt -j ${task.cpus} --times 2 -e 0.2 -a ${params.rev_primer} -o ${sampleid}_filtered2.fastq.gz ${sampleid}_filtered1.fastq.gz  >> ${sampleid}_cutadapt.log
    cutadapt -j ${task.cpus} --times 3 -e 0.2  -a "CCCGCGTACTCTGCGTTGATACCACTGCTT;min_overlap=14" -o ${sampleid}_filtered3.fastq.gz ${sampleid}_filtered2.fastq.gz  >> ${sampleid}_cutadapt.log
    cutadapt -j ${task.cpus} --times 2 -e 0.2 -g ${params.rev_primer_rc} -o ${sampleid}_cutadapt_filtered.fastq.gz ${sampleid}_filtered3.fastq.gz  >> ${sampleid}_cutadapt.log
    rm ${sampleid}_filtered*.fastq.gz

  elif [[ ${params.race3} == true ]]; then
    cutadapt -j ${task.cpus} --times 2 -e 0.2 -g ${params.fwd_primer} -o ${sampleid}_filtered1.fastq.gz ${sample} > ${sampleid}_cutadapt.log
    cutadapt -j ${task.cpus} --times 2 -e 0.2 -a ${params.fwd_primer_rc} -o ${sampleid}_cutadapt_filtered.fastq.gz ${sampleid}_filtered1.fastq.gz  >> ${sampleid}_cutadapt.log
    rm ${sampleid}_filtered*.fastq.gz
  fi
  """
}

process CHOPPER {
  publishDir "${params.outdir}/${sampleid}/canu", pattern:'*_filtered.fastq.gz', mode: 'link'
  publishDir "${params.outdir}/${sampleid}/canu", pattern: '*_chopper.log', mode: 'link'
  tag "${sampleid}"
  label 'large'

  container 'quay.io/biocontainers/chopper:0.5.0--hdcf5f25_2'

  input:
    tuple val(sampleid), path(sample), path(reference)

  output:
  path("${sampleid}_chopper.log")
    path("${sampleid}_filtered.fastq.gz")
    tuple val(sampleid), path("${sampleid}_filtered.fastq.gz"), path(reference), emit: chopper_filtered_fq

  script:
  """
  gunzip -c ${sample} | chopper -q ${params.chopper_qual_threshold} -l ${params.chopper_min_read_length} 2> ${sampleid}_chopper.log | gzip > ${sampleid}_filtered.fastq.gz
  """
}

process NANOFILT {
  publishDir "${params.outdir}/${sampleid}/canu", pattern:'*_filtered.fastq.gz', mode: 'link'
  //publishDir "${params.outdir}/${sampleid}/canu", pattern: '*_nanofilt.log', mode: 'link'
  tag "${sampleid}"
  label 'small'

  container 'quay.io/biocontainers/nanofilt:2.8.0--py_0'

  input:
    tuple val(sampleid), path(sample), path(reference)

  output:
    path("${sampleid}_filtered.fastq.gz")
    //path("${sampleid}_nanofilt.log")
    tuple val(sampleid), path("${sampleid}_filtered.fastq.gz"), path(reference), emit: nanofilt_filtered_fq

  script:
  """
  gunzip -c ${sample} | NanoFilt -q ${params.nanofilt_qual_threshold} -l ${params.nanofilt_min_read_length} | gzip > ${sampleid}_filtered.fastq.gz
  """
}

process CANU_RACE {
  publishDir "${params.outdir}/${sampleid}/canu", pattern:'*_unassembled.fasta', mode: 'link'
  publishDir "${params.outdir}/${sampleid}/canu", pattern:'*_assembly.fasta', mode: 'link'
  tag "${sampleid}"
  label 'xlarge'

  container 'quay.io/biocontainers/canu:2.2--ha47f30e_0'

  input:
    tuple val(sampleid), path(sample), path(reference)

  output:
    path("${sampleid}_canu_assembly.fasta")
    path("${sampleid}_canu_unassembled.fasta")
    tuple val(sampleid), path("${sampleid}_canu_assembly.fasta"), path(reference), emit: canu_race_assembly
    
  script:
  """
  canu -assemble -p ${sampleid} -d ${sampleid} \
    -corrected -trimmed \
    genomeSize=${params.canu_genome_size} \
    useGrid=false minOverlapLength=50  minReadLength=500 stopOnLowCoverage=0 corMinCoverage=0 \
    -nanopore ${sample} \
    

  cp ${sampleid}/${sampleid}.contigs.fasta ${sampleid}_canu_assembly.fasta
  cp ${sampleid}/${sampleid}.unassembled.fasta ${sampleid}_canu_unassembled.fasta
  """
}

process FLYE {
  publishDir "${params.outdir}/${sampleid}/flye", mode: 'link'
  tag "${sampleid}"
  label 'large'
  errorStrategy 'ignore'

  container "quay.io/biocontainers/flye:2.9.1--py310h590eda1_0"

  input:
    tuple val(sampleid), path(sample), path(reference)
  output:
    path 'outdir/*'
    path("${sampleid}_flye_assembly.fasta")
    tuple val(sampleid), path("${sampleid}_flye_assembly.fasta"), path(reference), emit: assembly
  script:
  """
  flye  --out-dir outdir --threads ${task.cpus} --read-error ${params.flye_read_error} --${params.flye_ont_mode} ${sample}
  
  if [[ ! -s outdir/assembly.fasta ]]
    then
        touch ${sampleid}_flye_assembly.fasta
  else 
    cp outdir/assembly.fasta ${sampleid}_flye_assembly.fasta
  fi
  """
}

process BLASTN {
  publishDir "${params.outdir}/${sampleid}/blastn", mode: 'link'
  tag "${sampleid}"
  label 'small'

  container 'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0'

  input:
    tuple val(sampleid), path(assembly), path(reference)
  output:
    path "BLASTN_${reference}_vs_${assembly}.txt"

  script:
  """
  blastn -query ${assembly} -subject ${reference} -evalue 1e-5 -out blastn_${reference}_vs_${assembly}.txt \
    -outfmt '6 qseqid sacc length pident mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp qcovs'
  
  echo "qseqid sacc length pident mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp qcovs" > header
  
  cat header blastn_${reference}_vs_${assembly}.txt > BLASTN_${reference}_vs_${assembly}.txt
  """
}

process MINIMAP2 {
  publishDir "${params.outdir}/${sampleid}/minimap2", mode: 'link'
  tag "${sampleid}"
  label 'medium'

  container 'quay.io/biocontainers/minimap2:2.24--h7132678_1'

  input:
    tuple val(sampleid), path(sample),path(reference)
  output:
    tuple val(sampleid), file("${sampleid}_aln.sam"), path(reference), emit: aligned_sample
  script:
  """
  minimap2 -a --MD ${reference} ${sample} > ${sampleid}_aln.sam
  """
}

process INFOSEQ {
  publishDir "${params.outdir}/${sampleid}/infoseq", mode: 'link'
  tag "${sampleid}"
  label 'small'

  container "quay.io/biocontainers/emboss:6.6.0--h1b6f16a_5"

  input:
    tuple val(sampleid), path(sample), path(reference)
  output:
    tuple val(sampleid), path(sample), path("${reference}_list.txt"), emit: infoseq_ref
  script:
  """
  infoseq ${reference} -only -name -length | sed 1d > ${reference}_list.txt
  """

}

process SAMTOOLS {
  publishDir "${params.outdir}/${sampleid}/samtools", mode: 'link'
  tag "${sampleid}"
  label 'small'

  container 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'

  input:
    tuple val(sampleid), path(sample), path(reference_list)
  output:
    tuple val(sampleid), path("${sampleid}_aln.sorted.bam"), path("${sampleid}_aln.sorted.bam.bai"), emit: sorted_sample
  script:
  """
  samtools view -bt ${reference_list} -o ${sampleid}_aln.bam ${sample}
  samtools sort -T /tmp/aln.sorted -o ${sampleid}_aln.sorted.bam ${sampleid}_aln.bam
  samtools index ${sampleid}_aln.sorted.bam
  """
}

process NANOQ {
  publishDir "${params.outdir}/${sampleid}/nano-q", mode: 'link'
  tag "${sampleid}"
  label 'medium'

  container 'ghcr.io/eresearchqut/nano-q:1.0.0'

  input:
    tuple val(sampleid), path(sorted_sample), path(sorted_sample_index)
  output:
    //path 'Results/*'
    path 'Results'

  script:
  """
  nano-q.py -b ${sorted_sample} -c ${params.nanoq_code_start} -l ${params.nanoq_read_length} -nr ${params.nanoq_num_ref} -q ${params.nanoq_qual_threshhold} -j ${params.nanoq_jump}
  """
}

workflow {
  if (params.samplesheet) {
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), file(row.sample_files), file(row.reference)) }
      .set{ ch_sample }
  } else { exit 1, "Input samplesheet file not specified!" }

  MERGE ( ch_sample )

  if ( params.canu)  {
    if ( params.race3 || params.race5 ) {
      /*
      if (params.cutadapt) {
        CUTADAPT_RACE ( MERGE.out.merged )
        }
      */
      CUTADAPT_RACE ( MERGE.out.merged )
      if (params.nanofilt) {
        NANOFILT ( CUTADAPT_RACE.out.cutadapt_filtered )
        CANU_RACE ( NANOFILT.out.nanofilt_filtered_fq )
        }
      else if (params.chopper) {
        CHOPPER ( CUTADAPT_RACE.out.cutadapt_filtered )
        CANU_RACE ( CHOPPER.out.chopper_filtered_fq )
        }
      BLASTN ( CANU_RACE.out.canu_race_assembly )
    }  
  }
  NANOPLOT ( ch_sample )
  if (params.flye) {
    FLYE ( MERGE.out.merged )
    BLASTN ( FLYE.out.assembly )
  }
  MINIMAP2 ( MERGE.out.merged )
  if (params.infoseq) {
    INFOSEQ ( MINIMAP2.out.aligned_sample )
    SAMTOOLS ( INFOSEQ.out.infoseq_ref )
    NANOQ ( SAMTOOLS.out.sorted_sample )
  }
}
