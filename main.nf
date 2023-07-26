#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage () {
    log.info """
    OViSP - ONT-based Viral Screening for Plants
    Marie-Emilie Gauthier 23/05/2023
    Craig Windell
    Roberto Barrero

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
        sampleid,sample_files
        SAMPLE01,/user/folder/sample.fastq.gz
        SAMPLE02,/user/folder/*.fastq.gz

        sample_files can refer to a folder with a number of
        files that will be merged in the pipeline

        --skip_porechop                 Skip porechop step
                              Default:  false
        --skip_nanofilt                 Skip nanofilt step
                              Default:  false
        --nanofilt_qual_threshold       Quality threshold
                              Default:  10
        --nanofilt_min_read_length      Length cut off for read size
                              Default:  500
        --skip_host_filtering           Skip host filtering step
                              Default:  false
        --plant_host_fasta              Fasta file of nucleotide sequences to filter
                                        null
        --skip_denovo_assembly          Skip de novo assembly step
                              Default:  false
        --canu                          Use Canu for de novo assembly step
                              Default:  false
        --canu_genome_size              Target genome size
                              Default:  '0.01m'
        --canu_options                  Canu options
                              Default:  ''
        --flye                          Use Flye for de novo assembly step
        --flye_ont_mode                 Select from nano-raw, nano-corr, nano-hq
                              Default:  'nano-raw'
        --flye_options                  Flye options
                              Default:  ''
        --skip_clustering               Skip clustering step using Rattle
                              Default:  false
        --rattle_min_len                Minimum length cut off for read size
                              Default:  250
        --rattle_max_len                Maximum length cut off for read size
                              Default:  2000

    """.stripIndent()
}
// Show help message
if (params.help) {
    helpMessage()
    exit 0
}
if (params.blastn_db != null) {
    blastn_db_name = file(params.blastn_db).name
    blastn_db_dir = file(params.blastn_db).parent
}
if (params.reference != null) {
    reference_name = file(params.reference).name
    reference_dir = file(params.reference).parent

}
if (params.kaiju_nodes != null & params.kaiju_dbname != null & params.kaiju_names != null) {
    kaiju_dbs_dir = file(params.kaiju_nodes).parent
}
if (params.krkdb != null) {
    krkdb_dir = file(params.krkdb).parent
}

switch (workflow.containerEngine) {
  case "singularity":
    bindbuild = "";
    if (params.blastn_db != null) {
      bindbuild = (bindbuild + "-B ${blastn_db_dir} ")
    }
    if (params.reference != null) {
      bindbuild = (bindbuild + "-B ${reference_dir} ")
    }
    if (params.kaiju_nodes != null & params.kaiju_dbname != null & params.kaiju_names != null) {
      bindbuild = (bindbuild + "-B ${kaiju_dbs_dir} ")
    }
    if (params.krkdb != null) {
      bindbuild = (bindbuild + "-B ${krkdb_dir} ")
    }
    bindOptions = bindbuild;
    break;
  default:
    bindOptions = "";
}

process MERGE {
  //publishDir "${params.outdir}/${sampleid}/merge", pattern: '*.fastq.gz', mode: 'link'
  tag "${sampleid}"
  label 'small'

  input:
    tuple val(sampleid), path(lanes)
  output:
    tuple val(sampleid), path("${sampleid}.fastq.gz"), emit: merged
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
    tuple val(sampleid), path(sample)
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
  publishDir "${params.outdir}/${sampleid}/cutadapt", pattern: '*_cutadapt.log', mode: 'link'
  tag "${sampleid}"
  label 'medium'

  container 'quay.io/biocontainers/cutadapt:4.1--py310h1425a21_1'

  input:
    tuple val(sampleid), path(sample)
  output:
    //path("${sampleid}_cutadapt_filtered.fastq.gz")
    path("${sampleid}_cutadapt.log")
    tuple val(sampleid), path("${sampleid}_cutadapt_filtered.fastq.gz"), emit: cutadapt_filtered
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
/*
process CHOPPER {
  //publishDir "${params.outdir}/${sampleid}/chopper", pattern:'*_filtered.fastq.gz', mode: 'link'
  publishDir "${params.outdir}/${sampleid}/chopper", pattern: '*_chopper.log', mode: 'link'
  tag "${sampleid}"
  label 'large'

  container 'quay.io/biocontainers/chopper:0.5.0--hdcf5f25_2'

  input:
    tuple val(sampleid), path(sample)

  output:
    path("${sampleid}_chopper.log")
    path("${sampleid}_filtered.fastq.gz")
    tuple val(sampleid), path("${sampleid}_filtered.fastq.gz"), emit: chopper_filtered_fq

  script:
  """
  gunzip -c ${sample} | chopper -q ${params.chopper_qual_threshold} -l ${params.chopper_min_read_length} 2> ${sampleid}_chopper.log | sed 's/ /_/g' | gzip > ${sampleid}_filtered.fastq.gz
  """
}
*/
process NANOFILT {
  publishDir "${params.outdir}/${sampleid}/nanofilt", pattern:'*_filtered.fastq.gz', mode: 'link'
  //publishDir "${params.outdir}/${sampleid}/canu", pattern: '*_nanofilt.log', mode: 'link'
  tag "${sampleid}"
  label 'small'

  container 'quay.io/biocontainers/nanofilt:2.8.0--py_0'

  input:
    tuple val(sampleid), path(sample)

  output:
    path("${sampleid}_filtered.fastq.gz")
    //path("${sampleid}_nanofilt.log")
    tuple val(sampleid), path("${sampleid}_filtered.fastq.gz"), emit: nanofilt_filtered_fq

  script:
  """
  gunzip -c ${sample} | NanoFilt -q ${params.nanofilt_qual_threshold} -l ${params.nanofilt_min_read_length} | gzip > ${sampleid}_filtered.fastq.gz
  """
}

process CANU {
  publishDir "${params.outdir}/${sampleid}/denovo", mode: 'link', overwrite: true
  tag "${sampleid}"
  memory "24GB"
  cpus "4"

  container 'quay.io/biocontainers/canu:2.2--ha47f30e_0'

  input:
    tuple val(sampleid), path(fastq)

  output:
    path("${sampleid}_canu_assembly.fasta")
    path("${sampleid}.canu.log")
    tuple val(sampleid), path("${sampleid}_canu.fastq"), path("${sampleid}_canu_assembly.fasta"), emit: assembly
    tuple val(sampleid), path("${sampleid}_canu_assembly.fasta"), emit: assembly2

    
  script:
  def canu_options = (params.canu_options) ? " ${params.canu_options}" : ''

  """
  canu -p ${sampleid} -d ${sampleid} \
    genomeSize=${params.canu_genome_size} \
    -nanopore ${fastq} ${canu_options} 2> ${sampleid}.canu.log

  if [[ ! -s ${sampleid}/${sampleid}.contigs.fasta ]]
    then
      touch ${sampleid}.canu_assembly.fasta
  else 
    cat ${sampleid}/${sampleid}.contigs.fasta ${sampleid}/${sampleid}.unassembled.fasta > ${sampleid}_canu_assembly.fasta
  fi
  cp ${fastq} ${sampleid}_canu.fastq
  """
}

process FLYE {
  publishDir "${params.outdir}/${sampleid}/denovo", mode: 'link'
  tag "${sampleid}"
  label 'large'
  

  container "quay.io/biocontainers/flye:2.9.1--py310h590eda1_0"

  input:
    tuple val(sampleid), path(sample)
  output:
    path 'outdir/*'
    path("${sampleid}_flye_assembly.fasta")
    tuple val(sampleid), path("${sampleid}_flye_assembly.fasta"), emit: assembly
  script:
  
  """
  flye  --out-dir outdir --threads ${task.cpus}  --${params.flye_ont_mode} ${sample}
  
  if [[ ! -s outdir/assembly.fasta ]]
    then
        touch ${sampleid}_flye_assembly.fasta
  else 
    cp outdir/assembly.fasta ${sampleid}_flye_assembly.fasta
  fi
  """
}
/*
errorStrategy 'ignore'
flye  --out-dir outdir --threads ${task.cpus} --read-error ${params.flye_read_error} --${params.flye_ont_mode} ${sample}

*/

process BLASTN2REF {
  publishDir "${params.outdir}/${sampleid}/blast_to_ref", mode: 'link'
  tag "${sampleid}"
  label 'small'
  containerOptions "${bindOptions}"

  container 'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0'

  input:
    tuple val(sampleid), path(assembly)
  output:
    path "BLASTN_reference_vs_${assembly}.txt"

  script:
  """
  blastn -query ${assembly} -subject ${reference_dir}/${reference_name} -evalue 1e-3 -out blastn_reference_vs_${assembly}.txt \
  -outfmt '6 qseqid sacc length pident mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp qcovs' -max_target_seqs 5

  echo "qseqid sacc length pident mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp qcovs" > header
  cat header blastn_reference_vs_${assembly}.txt > BLASTN_reference_vs_${assembly}.txt

  """
}
/*
process MINIMAP2 {
  publishDir "${params.outdir}/${sampleid}/denovo", mode: 'link'
  tag "${sampleid}"
  label 'large'
  containerOptions "${bindOptions}"

  container 'quay.io/biocontainers/minimap2:2.24--h7132678_1'

  input:
    tuple val(sampleid), path(fastq)
  output:
    tuple val(sampleid), path(fastq), path("${sampleid}_minimap.paf"), emit: paf
  script:
  """
  minimap2 -x ava-ont -t ${params.minimap_threads} ${fastq} ${fastq}  > ${sampleid}_minimap.paf
  """
}

process MINIASM {
  publishDir "${params.outdir}/${sampleid}/denovo", mode: 'link'
  tag "${sampleid}"
  label 'large'
  containerOptions "${bindOptions}"

  container 'quay.io/biocontainers/miniasm:0.3--he4a0461_2'

  input:
    tuple val(sampleid), path(fastq), path(paf)
  output:
    file("${sampleid}_miniasm.fasta")
  script:
  """
  miniasm -f ${fastq} ${paf} > ${sampleid}_miniasm.gfa
  awk '/^S/{print ">"\$2"\\n"\$3}' ${sampleid}_miniasm.gfa > ${sampleid}_miniasm.fasta
  """
}
*/
process MINIMAP2_REF {
  publishDir "${params.outdir}/${sampleid}/minimap2", mode: 'link'
  tag "${sampleid}"
  label 'medium'
  containerOptions "${bindOptions}"

  container 'quay.io/biocontainers/minimap2:2.24--h7132678_1'

  input:
    tuple val(sampleid), path(sample)
  output:
    tuple val(sampleid), file("${sampleid}_aln.sam"), emit: aligned_sample
  script:
  """
  minimap2 -a --MD ${reference_dir}/${reference_name} ${sample} > ${sampleid}_aln.sam
  """
}

process INFOSEQ {
  publishDir "${params.outdir}/${sampleid}/infoseq", mode: 'link'
  tag "${sampleid}"
  label 'small'
  containerOptions "${bindOptions}"

  container "quay.io/biocontainers/emboss:6.6.0--h1b6f16a_5"

  input:
    tuple val(sampleid), path(sample)
  output:
    tuple val(sampleid), path(sample), emit: infoseq_ref
  script:
  """
  infoseq ${reference_dir}/${reference_name} -only -name -length | sed 1d > ${reference_name}_list.txt
  """
}

process SAMTOOLS {
  publishDir "${params.outdir}/${sampleid}/samtools", mode: 'link'
  tag "${sampleid}"
  label 'small'

  container 'quay.io/biocontainers/samtools:1.16.1--h6899075_1'

  input:
    tuple val(sampleid), path(sample)
  output:
    tuple val(sampleid), path("${sampleid}_aln.sorted.bam"), path("${sampleid}_aln.sorted.bam.bai"), emit: sorted_sample
  script:
  """
  samtools view -bt ${reference_dir}/${reference_name} -o ${sampleid}_aln.bam ${sample}
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
    tuple val(sampleid), path(sorted_sample)
  output:
    //path 'Results/*'
    path 'Results'

  script:
  """
  nano-q.py -b ${sorted_sample} -c ${params.nanoq_code_start} -l ${params.nanoq_read_length} -nr ${params.nanoq_num_ref} -q ${params.nanoq_qual_threshhold} -j ${params.nanoq_jump}
  """
}
/*
process PORECHOP {
	tag "${sampleid}"
	label "xlarge2"
	publishDir "$params.outdir/${sampleid}/porechop",  mode: 'copy', pattern: "*.log", saveAs: { filename -> "${sample}_$filename" }

  container = 'docker://quay.io/biocontainers/porechop:0.2.3_seqan2.1.1--0'

	input:
		tuple val(sampleid), path(sample)
	output:
		tuple val(sampleid), file("porechop_trimmed.fastq.gz"), emit: porechop_trimmed_fq
	script:
	"""
	porechop -i ${sample} -t ${params.porechop_threads} -o porechop_trimmed.fastq.gz ${params.porechop_args}
	"""
}
*/
process PORECHOP_ABI {
  tag "${sampleid}"
  label "xlarge2"
  publishDir "$params.outdir/${sampleid}/porechop",  mode: 'copy'

  container = 'docker://quay.io/biocontainers/porechop_abi:0.5.0--py38he0f268d_2'

  input:
    tuple val(sampleid), path(sample)

  output:
    tuple val(sampleid), file("${sampleid}_porechop_trimmed.fastq.gz"), emit: porechopabi_trimmed_fq

  script:
  """
  porechop_abi -abi -i ${sample} -t ${params.porechop_threads} -o ${sampleid}_porechop_trimmed.fastq.gz ${params.porechop_args}
  """
}

process REFORMAT {
  tag "${sampleid}"
  label "large"
  publishDir "$params.outdir/${sampleid}",  mode: 'copy'

  container = 'docker://quay.io/biocontainers/bbmap:39.01--h92535d8_1'

  input:
  tuple val(sampleid), path(fastq)
  output:
  tuple val(sampleid), path("${sampleid}_quality_trimmed.fastq.gz"), emit: reformatted_fq

  script:
  """
  reformat.sh in=${fastq} out=${sampleid}_quality_trimmed.fastq.gz trd
  """
}

process CAP3 {
  tag "${sampleid}"
  label "large"
  time "3h"
  publishDir "$params.outdir/$sampleid/cap3", mode: 'copy', pattern: '*_cap3.fasta', saveAs: { filename -> "${sampleid}_cap3.fasta"}

  container = 'docker://quay.io/biocontainers/cap3:10.2011--h779adbc_3'

  input:
  tuple val(sampleid), path(fasta)
  output:
  tuple val(sampleid), path("${sampleid}_cap3.fasta"), emit: contigs

  script:
  """
  cap3 ${fasta}
  cat ${fasta}.cap.singlets ${fasta}.cap.contigs > ${sampleid}_cap3.fasta
  """
}
/*
process BLASTN {
  cpus "${params.blast_threads}"
  tag "${sampleid}"
  label "xlarge"
  time "5h"
  containerOptions "${bindOptions}"
  publishDir "$params.outdir/$sampleid/blast",  mode: 'link', overwrite: true, pattern: '*.bls', saveAs: { filename -> "${sampleid}_blastn_vs_NT.bls"}

  container 'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0'

  input:
  tuple val(sampleid), path(assembly)
  output:
  path("*.bls")
  tuple val(sampleid), path("${sampleid}_${params.blastn_method}_vs_NT.bls"), emit: blast_results

  script:
  def blast_task_param = (params.blastn_method == "blastn") ? "-task blastn" : ''
  """
  cp ${blastn_db_dir}/taxdb.btd .
  cp ${blastn_db_dir}/taxdb.bti .
  blastn ${blast_task_param} \
    -query ${assembly} \
    -db ${params.blastn_db} \
    -out ${sampleid}_${params.blastn_method}_vs_NT.bls \
    -evalue 1e-3 \
    -num_threads ${params.blast_threads} \
    -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames' \
    -max_target_seqs 5
"""
}
*/
process EXTRACT_VIRAL_BLAST_HITS {
  tag "${sampleid}"
  label "medium"
  publishDir "$params.outdir/$sampleid/blastn",  mode: 'link', overwrite: true

  container = 'docker://infrahelpers/python-light:py310-bullseye'

  input:
  tuple val(sampleid), path(blast_results)
  output:
  file "${sampleid}_megablast_vs_NT_top_hits.txt"
  file "${sampleid}_megablast_vs_NT_top_viral_hits.txt"
  file "${sampleid}_megablast_vs_NT_top_viral_spp_hits.txt"

  script:
  """  
  cat ${blast_results} > ${sampleid}_megablast_vs_NT.txt

  select_top_blast_hit.py --sample_name ${sampleid} --megablast_results ${sampleid}_megablast_vs_NT.txt
  """
}

process CONCATENATE_FASTA {
  tag "${sampleid}"
  label "medium"
  publishDir "${params.outdir}/${sampleid}", mode: 'link'

  container = 'docker://quay.io/biocontainers/seqtk:1.3--h5bf99c6_3'

  input:
  tuple val(sampleid), path("${sampleid}_canu_assembly.fasta")
  tuple val(sampleid), path("${sampleid}_cap3.fasta")
  tuple val(sampleid), path("${sampleid}.fasta")
  output:
  file "${sampleid}_merged.fasta"
  tuple val(sampleid), path("*_merged.fasta"), emit: assembly

  script:
  """
  seqtk seq -l0 ${sampleid}_canu_assembly.fasta > ${sampleid}_canu_assembly_1l.fasta
  seqtk seq -l0 ${sampleid}_cap3.fasta >  ${sampleid}_cap3_1l.fasta
  seqtk seq -l0 ${sampleid}.fasta >  ${sampleid}_1l.fasta 
  cat  ${sampleid}_canu_assembly_1l.fasta ${sampleid}_cap3_1l.fasta  ${sampleid}.fasta > ${sampleid}_merged.fasta
  """

}

process BLASTN_SPLIT {
  publishDir "${params.outdir}/${sampleid}/blastn", mode: 'link'
  tag "${sampleid}"
  containerOptions "${bindOptions}"
  time "12h"
  memory "48GB"
  cpus "4"

  container 'quay.io/biocontainers/blast:2.13.0--hf3cf87c_0'

  input:
    tuple val(sampleid), path(assembly)
  output:
    tuple val(sampleid), path("${sampleid}*_blastn_vs_NT.bls"), emit: blast_results

  script:
  def blastoutput = assembly.getBaseName() + "_blastn_vs_NT.bls"
  """
  cp ${blastn_db_dir}/taxdb.btd .
  cp ${blastn_db_dir}/taxdb.bti .
  blastn -query ${assembly} \
    -db ${params.blastn_db} \
    -out ${blastoutput} \
    -evalue 1e-3 \
    -num_threads ${params.blast_threads} \
    -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames' \
    -max_target_seqs 5
  """
}
/*
KAIJU notes
The default run mode is Greedy with three allowed mismatches. 
The number of allowed mismatches can be changed using option -e.
In Greedy mode, matches are filtered by a minimum length and score and their E-value (similar to blastp)
The cutoffs for minimum required match length and match score can be changed using the options -m (default: 11) and -s (default: 65)
Minimum E-value can be adjusted with the option -E (default 0.01).

For fastest classification, , use MEM mode with option '-a mem' and multiple parallel threads (-z)
Greedy run mode yields a higher sensitivity compared with MEM mode.
For lowest memory usage use the proGenomes reference database. The number of parallel threads has only little impact on memory usage.


Further, the choice of the minimum required match length (-m) in MEM mode or match score (-s) in Greedy mode governs the trade-off between 
sensitivity and precision of the classification. Please refer to the paper for a discussion on this topic.

Option -x enables filtering of query sequences containing low-complexity regions by using the SEG algorithm from the blast+ package. 
It is enabled by default and can be disabled by the -X option. SEG filtering is always recommended in order to avoid 
false positive taxon assignments that are caused by spurious matches due to simple repeat patterns or other sequencing noise.

The accuracy of the classification depends both on the choice of the reference database and the chosen options when running Kaiju. 
These choices also affect the speed and memory usage of Kaiju.

For highest sensitivity, it is recommended to use the nr database (+eukaryotes) as a reference database because it is the most comprehensive 
set of protein sequences. Alternatively, use proGenomes over Refseq for increased sensitivity.
NOTE, viroid will not be detected using this approach

Mandatory arguments:

  -t FILENAME   Name of nodes.dmp file
  -f FILENAME   Name of database (.fmi) file
  -i FILENAME   Name of input file containing reads in FASTA or FASTQ format
Optional arguments:

  -j FILENAME   Name of second input file for paired-end reads
  -o FILENAME   Name of output file. If not specified, output will be printed to STDOUT
  -z INT        Number of parallel threads for classification (default: 1)
  -a STRING     Run mode, either "mem"  or "greedy" (default: greedy)
  -e INT        Number of mismatches allowed in Greedy mode (default: 3) #set to 1 in kodoja
  -m INT        Minimum match length (default: 11) #set to 15 in kodoja
  -s INT        Minimum match score in Greedy mode (default: 65) #set to 85 in Kodoja
  -E FLOAT      Minimum E-value in Greedy mode
  -x            Enable SEG low complexity filter (enabled by default)
  -X            Disable SEG low complexity filter
  -p            Input sequences are protein sequences
  -v            Enable verbose output
*/
process KAIJU {
  publishDir "${params.outdir}/${sampleid}/kaiju", mode: 'link'
  label 'process_high'
  container 'quay.io/biocontainers/kaiju:1.8.2--h5b5514e_1'
  containerOptions "${bindOptions}"

  input:
  tuple val(sampleid), path(fastq)

  output:
  tuple val(sampleid), path('*kaiju.krona'), emit: results
  file "${sampleid}_kaiju_name.tsv"
  file "${sampleid}_kaiju_summary.tsv"
  file "${sampleid}_kaiju.krona"

  script:
  """

  kaiju \\
      -z ${params.kaiju_threads} \\
      -t ${params.kaiju_nodes}  \\
      -f ${params.kaiju_dbname} \\
      -o ${sampleid}_kaiju.tsv \\
      -i ${fastq} \\
      -v
  
  kaiju-addTaxonNames -t ${params.kaiju_nodes} -n ${params.kaiju_names} -i ${sampleid}_kaiju.tsv -o ${sampleid}_kaiju_name.tsv
  kaiju2table -t ${params.kaiju_nodes} -r species -n ${params.kaiju_names} -o ${sampleid}_kaiju_summary.tsv ${sampleid}_kaiju.tsv
  kaiju2krona -t ${params.kaiju_nodes} -n ${params.kaiju_names} -i ${sampleid}_kaiju.tsv -o ${sampleid}_kaiju.krona
  
  """
}

process KRONA {
  publishDir "${params.outdir}/${sampleid}/krona", mode: 'link'
  label 'large'
  container 'quay.io/biocontainers/krona:2.8.1--pl5321hdfd78af_1'
  containerOptions "${bindOptions}"

  input:
  tuple val(sampleid), path(krona_input)

  output:
  file "${sampleid}_krona.html"

  script:
  """
  ktImportText -o ${sampleid}_krona.html ${krona_input}
  """
}

/*

KRAKEN2
https://github.com/DerrickWood/kraken2/blob/master/docs/MANUAL.markdown
If you get the Loading database information...classify: Error reading in hash table error allocate more memory to the job. Keep incrementing the memory request until you see a status message in the job log like

Usage: kraken2 [options] <filename(s)>

  Options:
    --db NAME               Name for Kraken 2 DB
                            (default: none)
    --threads NUM           Number of threads (default: 1)
    --quick                 Quick operation (use first hit or hits)
    --unclassified-out FILENAME
                            Print unclassified sequences to filename
    --classified-out FILENAME
                            Print classified sequences to filename
    --output FILENAME       Print output to filename (default: stdout); "-" will
                            suppress normal output
    --confidence FLOAT      Confidence score threshold (default: 0.0); must be
                            in [0, 1].
    --minimum-base-quality NUM
                            Minimum base quality used in classification (def: 0,
                            only effective with FASTQ input).
    --report FILENAME       Print a report with aggregrate counts/clade to file
    --use-mpa-style         With --report, format report output like Kraken 1's
                            kraken-mpa-report
    --report-zero-counts    With --report, report counts for ALL taxa, even if
                            counts are zero
    --report-minimizer-data With --report, report minimizer and distinct minimizer
                            count information in addition to normal Kraken report
    --memory-mapping        Avoids loading database into RAM
    --paired                The filenames provided have paired-end reads
    --use-names             Print scientific names instead of just taxids
    --gzip-compressed       Input files are compressed with gzip
    --bzip2-compressed      Input files are compressed with bzip2
    --minimum-hit-groups NUM
                            Minimum number of hit groups (overlapping k-mers
                            sharing the same minimizer) needed to make a call
                            (default: 3)
    --help                  Print this message
    --version               Print version information

32 min using 4 cpus and 164 Gb of mem, test with 2 cpus instead
*/

process KRAKEN2 {
	tag "${sampleid}"
	label 'process_high2'
	publishDir "$params.outdir/$sampleid/kraken",  mode: 'link'
  container 'quay.io/biocontainers/kraken2:2.1.3--pl5321hdcf5f25_0'
  containerOptions "${bindOptions}"

	input:
		tuple val(sampleid), path(fastq)

	output:
		file("${sampleid}.kraken2")
    file("${sampleid}_kraken_report.txt")
    file("${sampleid}_seq_ids.txt")
    tuple val(sampleid), path("${sampleid}_kraken_report.txt"), emit: results
	script:
	"""
	kraken2 --db ${params.krkdb} \\
          --use-names \\
          --gzip-compressed \\
          --threads 2 \\
          --report ${sampleid}_kraken_report.txt  \\
          --report-minimizer-data \\
          --minimum-hit-groups 3 \\
          ${fastq} > ${sampleid}.kraken2

  #extract the reads IDs
  echo "seq_id" > ${sampleid}_seq_ids.txt
	awk -F "\\t" '{print \$2}' ${sampleid}.kraken2 >> ${sampleid}_seq_ids.txt
  """
}
/*
gawk 'match($0, pattern, ary) {print ary[1]}'
	echo "seq_id" > seq_ids.txt 
	awk -F "\\t" '{print \$2}' krakenreport.txt >> seq_ids.txt       
	gawk -F "\\t" 'match(\$0, /\\(taxid\s([0-9]+)\\)/, ary) {print ary[1]}' krakenreport.txt | taxonkit lineage --data-dir $taxondb > lineage.txt
	cat lineage.txt | taxonkit reformat --data-dir $taxondb | csvtk -H -t cut -f 1,3 | csvtk -H -t sep -f 2 -s ';' -R > seq_tax.txt
	cat lineage.txt | taxonkit reformat -P --data-dir $taxondb | csvtk -H -t cut -f 1,3 > seq_tax_otu.txt
	paste seq_ids.txt seq_tax.txt > kraken_report_annotated.txt
	paste seq_ids.txt seq_tax_otu.txt > kraken_report_annotated_otu.txt
	"""

set +eu
	sed '/^@/s/.\s./_/g' ${filtered} > krkinput.fastq
*/

process BRACKEN {
  tag "${sampleid}"
	label 'medium'
	publishDir "$params.outdir/$sampleid/bracken",  mode: 'link'
  container 'quay.io/biocontainers/bracken:2.8--py310h0dbaff4_1'
  containerOptions "${bindOptions}"

	input:
		tuple val(sampleid), path(kraken_report)

	output:
		file("${sampleid}_bracken_report.txt")
	script:
	"""
	est_abundance.py -i ${kraken_report} \\
                  -k ${params.krkdb}/database50mers.kmer_distrib \\
                  -t 1 \\
                  -l S -o ${sampleid}_bracken_report.txt
  """
}


process RATTLE {
  publishDir "${params.outdir}/${sampleid}/clustering", mode: 'link'
  tag "${sampleid}"
  label 'medium'
  containerOptions "${bindOptions}"

  input:
  tuple val(sampleid), path(fastq)

  output:
  tuple val(sampleid), path("transcriptome.fq"), emit: clusters
  tuple val(sampleid), path(fastq), path("transcriptome.fq"), emit: clusters2

  container =  'ghcr.io/eresearchqut/rattle-image:0.0.1'

  script:
  """
  rattle cluster -i ${fastq} -t 2 --lower-length ${params.rattle_min_len} --upper-length ${params.rattle_max_len}  -o . --rna
  rattle cluster_summary -i ${fastq} -c clusters.out > cluster_summary.txt
  mkdir clusters
  rattle extract_clusters -i ${fastq} -c clusters.out --fastq -o clusters
  rattle correct -i ${fastq} -c clusters.out -t 2
  rattle polish -i consensi.fq -t 2 --rna --summary
  """

}

include { MINIMAP2_ALIGN as FILTER_HOST} from './modules.nf'
include { EXTRACT_READS as EXTRACT_READS_STEP1 } from './modules.nf'
include { EXTRACT_READS as EXTRACT_READS_STEP2 } from './modules.nf'
include { EXTRACT_READS as EXTRACT_READS_STEP3 } from './modules.nf'
include { MINIMAP2_ALIGN as MAP_BACK_TO_ASSEMBLY} from './modules.nf'
include { MINIMAP2_ALIGN as MAP_BACK_TO_CLUSTERS } from './modules.nf'
include { FASTQ2FASTA as FASTQ2FASTA_STEP1} from './modules.nf'
include { FASTQ2FASTA as FASTQ2FASTA_STEP2} from './modules.nf'

workflow {
  if (params.samplesheet) {
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), file(row.sample_files)) }
      .set{ ch_sample }
  } else { exit 1, "Input samplesheet file not specified!" }
  
  MERGE ( ch_sample )
  NANOPLOT ( MERGE.out.merged )

  // Data pre-processing
  if (!params.skip_porechop & !params.skip_nanofilt) {
    PORECHOP_ABI ( MERGE.out.merged )
    NANOFILT ( PORECHOP_ABI.out.porechopabi_trimmed_fq )
    REFORMAT( NANOFILT.out.nanofilt_filtered_fq )
  }
  else if (!params.skip_porechop & params.skip_nanofilt) {
    PORECHOP_ABI ( MERGE.out.merged )
    REFORMAT( PORECHOP_ABI.out.porechopabi_trimmed_fq )
  }
  else if (params.skip_porechop & !params.skip_nanofilt) {
    NANOFILT ( MERGE.out.merged )
    REFORMAT( NANOFILT.out.nanofilt_filtered_fq )
  }

  else if (params.skip_porechop & params.skip_nanofilt ) {
    REFORMAT(  MERGE.out.merged )
  }

  if (!params.skip_host_filtering & !params.skip_denovo_assembly & !params.skip_clustering) {
    FILTER_HOST( REFORMAT.out.reformatted_fq, params.plant_host_fasta )
    EXTRACT_READS_STEP1( FILTER_HOST.out.sequencing_ids )

    if (params.canu) {
      CANU( EXTRACT_READS_STEP1.out.unaligned_fq )
      MAP_BACK_TO_ASSEMBLY ( CANU.out.assembly )
      EXTRACT_READS_STEP2( MAP_BACK_TO_ASSEMBLY.out.sequencing_ids )
      RATTLE ( EXTRACT_READS_STEP2.out.unaligned_fq )
      FASTQ2FASTA_STEP1( RATTLE.out.clusters )
      CAP3( FASTQ2FASTA_STEP1.out.fasta )
      MAP_BACK_TO_CLUSTERS ( RATTLE.out.clusters2 )
      EXTRACT_READS_STEP3 ( MAP_BACK_TO_CLUSTERS.out.sequencing_ids )
      FASTQ2FASTA_STEP2( EXTRACT_READS_STEP3.out.unaligned_fq )
      CONCATENATE_FASTA(CANU.out.assembly2, CAP3.out.contigs, FASTQ2FASTA_STEP2.out.fasta)
      BLASTN_SPLIT( CONCATENATE_FASTA.out.assembly).splitFasta(by: 25000, file: true)
      BLASTN_SPLIT.out.blast_results
        .groupTuple()
        .set { ch_blastresults }
      EXTRACT_VIRAL_BLAST_HITS( ch_blastresults )
    }
    else if (params.flye) {
      FLYE( EXTRACT_READS_STEP1.out.unaligned_fq )
      MAP_BACK_TO_ASSEMBLY ( FLYE.out.assembly )
      EXTRACT_READS_STEP2( MAP_BACK_TO_ASSEMBLY.out.sequencing_ids )
      FASTQ2FASTA_STEP1( RATTLE.out.clusters )
      CAP3( FASTQ2FASTA_STEP1.out.fasta )
      BLASTN_SPLIT( FLYE.out.assembly.mix(CAP3.out.contigs).collect().splitFasta(by: 25000, file: true) )
    }
  }

  else if (params.skip_host_filtering & !params.skip_denovo_assembly & !params.skip_clustering) {
    if (params.canu) {
      CANU( REFORMAT.out.reformatted_fq )
      MAP_BACK_TO_ASSEMBLY ( CANU.out.assembly )
      EXTRACT_READS_STEP2( MAP_BACK_TO_ASSEMBLY.out.sequencing_ids )
      RATTLE ( EXTRACT_READS_STEP2.out.unaligned_fq )
      FASTQ2FASTA_STEP1( RATTLE.out.clusters )
      CAP3( FASTQ2FASTA_STEP1.out.fasta )
      BLASTN_SPLIT( CANU.out.assembly2.mix(CAP3.out.contigs).collect().splitFasta(by: 25000, file: true) )
    }
    else if (params.flye) {
      FLYE( REFORMAT.out.reformatted_fq )
      MAP_BACK_TO_ASSEMBLY ( FLYE.out.assembly )
      EXTRACT_READS_STEP2( MAP_BACK_TO_ASSEMBLY.out.sequencing_ids )
      RATTLE ( EXTRACT_READS_STEP2.out.unaligned_fq )
      FASTQ2FASTA_STEP1( RATTLE.out.clusters )
      CAP3( FASTQ2FASTA_STEP1.out.fasta )
      BLASTN_SPLIT( FLYE.out.assembly2.mix(CAP3.out.contigs).collect().splitFasta(by: 25000, file: true) )
    }
  }
  
  else if (params.skip_host_filtering & params.skip_denovo_assembly & !params.skip_clustering) {
    RATTLE ( REFORMAT.out.reformatted_fq )
    FASTQ2FASTA_STEP1( RATTLE.out.clusters )
    CAP3( FASTQ2FASTA_STEP1.out.fasta )
    
    if (params.blast_vs_ref) {
      BLASTN2REF ( CAP3.out.contigs )
    }
    else {
      BLASTN_SPLIT( CAP3.out.contigs.splitFasta(by: 25000, file: true) )
    }
  }

  else if (!params.skip_host_filtering & params.skip_denovo_assembly & !params.skip_clustering) {
    FILTER_HOST( REFORMAT.out.reformatted_fq, params.plant_host_fasta )
    EXTRACT_READS_STEP1( FILTER_HOST.out.sequencing_ids )

    RATTLE( EXTRACT_READS_STEP1.out.unaligned_fq )
    CAP3( RATTLE.out.clusters )
    
    if (params.blast_vs_ref) {
      BLASTN2REF ( CAP3.out.contigs )
    }
    else {
      BLASTN_SPLIT( CAP3.out.contigs.splitFasta(by: 25000, file: true) )
    }
  }

  else if (params.skip_host_filtering & !params.skip_denovo_assembly & params.skip_clustering) {
    if (params.canu) {
      CANU( REFORMAT.out.reformatted_fq )
      BLASTN_SPLIT( CANU.out.assembly2.splitFasta(by: 25000, file: true) )
      BLASTN_SPLIT.out.blast_results
        .groupTuple()
        .set { ch_blastresults }
      EXTRACT_VIRAL_BLAST_HITS( ch_blastresults )
    }
    else if (params.flye) {
      FLYE( REFORMAT.out.reformatted_fq )
      BLASTN_SPLIT( FLYE.out.assembly.splitFasta(by: 25000, file: true) )
      BLASTN_SPLIT.out.blast_results
        .groupTuple()
        .set { ch_blastresults }
      EXTRACT_VIRAL_BLAST_HITS( ch_blastresults )
    }
  }

  else if (!params.skip_host_filtering & !params.skip_denovo_assembly & params.skip_clustering) {
    FILTER_HOST( REFORMAT.out.reformatted_fq, params.plant_host_fasta )
    EXTRACT_READS_STEP1( FILTER_HOST.out.sequencing_ids )

    if (params.canu) {
      CANU( EXTRACT_READS_STEP1.out.unaligned_fq )
      BLASTN_SPLIT( CANU.out.assembly2.splitFasta(by: 25000, file: true) )
      BLASTN_SPLIT.out.blast_results
        .groupTuple()
        .set { ch_blastresults }
      EXTRACT_VIRAL_BLAST_HITS( ch_blastresults )
    }
    else if (params.flye) {
      FLYE( EXTRACT_READS_STEP1.out.unaligned_fq )
      BLASTN_SPLIT( FLYE.out.assembly.splitFasta(by: 25000, file: true) )
      BLASTN_SPLIT.out.blast_results
        .groupTuple()
        .set { ch_blastresults }
      EXTRACT_VIRAL_BLAST_HITS( ch_blastresults )
    }
  }
  
  //else if (params.skip_host_filtering & params.skip_denovo_assembly & params.skip_clustering) {
    //just perform direct alignment 
  if (params.map2ref) {
    MINIMAP2_REF ( REFORMAT.out.reformatted_fq )
    if (params.infoseq) {
      INFOSEQ ( MINIMAP2_REF.out.aligned_sample )
      SAMTOOLS ( INFOSEQ.out.infoseq_ref )
      NANOQ ( SAMTOOLS.out.sorted_sample )
    }
  }
  
   
  if (!params.skip_read_classification) {
    //just perform direct read search
    if (params.megablast) {
    FASTQ2FASTA_STEP1( REFORMAT.out.reformatted_fq )
    BLASTN_SPLIT( FASTQ2FASTA_STEP1.out.fasta.splitFasta(by: params.blast_split_factor, file: true) )
    BLASTN_SPLIT.out.blast_results
      .groupTuple()
      .set { ch_blastresults }
    EXTRACT_VIRAL_BLAST_HITS( ch_blastresults )
    }
    if (params.kaiju) {
      KAIJU ( REFORMAT.out.reformatted_fq )
      KRONA ( KAIJU.out.results)
    }
    if (params.kraken2) {
      KRAKEN2 ( REFORMAT.out.reformatted_fq )
      BRACKEN ( KRAKEN2.out.results )
    }
  }
}