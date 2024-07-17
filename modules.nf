if (params.blastn_db != null) {
    blastn_db_name = file(params.blastn_db).name
    blastn_db_dir = file(params.blastn_db).parent
}
if (params.host_fasta != null) {
    host_fasta_dir = file(params.host_fasta).parent
}

switch (workflow.containerEngine) {
  case "singularity":
    bindbuild = "";
    if (params.blastn_db != null) {
      bindbuild = (bindbuild + "-B ${blastn_db_dir} ")
    }
    if (params.host_fasta != null) {
      bindbuild = (bindbuild + "-B ${host_fasta_dir} ")
    }
    bindOptions = bindbuild;
    break;
  default:
    bindOptions = "";
}

process MINIMAP2_ALIGN_DNA {
  tag "${sampleid}"
  label "setting_8"

  input:
  tuple val(sampleid), path(fastq)
  path(reference)
  output:
  tuple val(sampleid), path(fastq), path("${sampleid}_unaligned_ids.txt"), emit: sequencing_ids

  script:
  """
  minimap2 -ax map-ont -L ${reference} ${fastq} -t ${task.cpus} > ${sampleid}.sam
  awk '\$6 == "*" { print \$0 }' ${sampleid}.sam | cut -f1 | uniq >  ${sampleid}_unaligned_ids.txt
  """
}

process MINIMAP2_ALIGN_RNA {
  tag "${sampleid}"
  label "setting_8"
  containerOptions "${bindOptions}"

  input:
  tuple val(sampleid), path(fastq)
  path(reference)
  output:
  tuple val(sampleid), path(fastq), path("${sampleid}_unaligned_ids.txt"), emit: sequencing_ids

  script:
  """
  minimap2 -ax splice -uf -k14 -L ${reference} ${fastq} -t ${task.cpus} > ${sampleid}.sam
  awk '\$6 == "*" { print \$0 }' ${sampleid}.sam | cut -f1 | sort | uniq >  ${sampleid}_unaligned_ids.txt
  """
}

process EXTRACT_READS {
  tag "${sampleid}"
  label "setting_11"
  publishDir "${params.outdir}/${sampleid}/host_filtering", mode: 'copy', pattern: '{*.fastq.gz,*reads_count.txt}'

  input:
  tuple val(sampleid), path(fastq), path(unaligned_ids)
  output:
  path("*reads_count.txt"), emit: read_counts
  file("${sampleid}_unaligned_reads_count.txt")
  file("${sampleid}_unaligned.fastq.gz")
  tuple val(sampleid), path("*_unaligned.fastq"), emit: unaligned_fq

  script:
  """
  seqtk subseq ${fastq} ${unaligned_ids} > ${sampleid}_unaligned.fastq
  gzip -c ${sampleid}_unaligned.fastq > ${sampleid}_unaligned.fastq.gz
  
  n_lines=\$(expr \$(cat ${sampleid}_unaligned.fastq | wc -l) / 4)
  echo \$n_lines > ${sampleid}_unaligned_reads_count.txt
  """
}


process FASTQ2FASTA {
  tag "${sampleid}"
  label "setting_2"

  input:
  tuple val(sampleid), path(fastq)
  output:
  tuple val(sampleid), path("${sampleid}.fasta"), emit: fasta

  script:
  """
  seqtk seq -A -C ${fastq} > ${sampleid}.fasta
  """
}

process NANOPLOT {
  publishDir "${params.outdir}/${sampleid}/qc/nanoplot",  pattern: '{*NanoPlot-report.html}', mode: 'link'
  publishDir "${params.outdir}/${sampleid}/qc/nanoplot",  pattern: '{*NanoStats.txt}', mode: 'link'
  publishDir "${params.outdir}/${sampleid}/qc/nanoplot",  pattern: '{*LengthvsQualityScatterPlot_dot.html}', mode: 'link'
  tag "${sampleid}"
  label "setting_10"

  input:
    tuple val(sampleid), path(sample)
  output:
    path("*NanoPlot-report.html")
    path("*NanoStats.txt")
    path("*LengthvsQualityScatterPlot_dot.html")
    path("*NanoStats.txt"), emit: read_counts
  
  script:
  """
  if [[ ${sample} == *trimmed.fastq.gz ]] || [[ ${sample} == *filtered.fastq.gz ]] ;
  then
    NanoPlot -t 8 --fastq ${sample} --prefix ${sampleid}_filtered_ --plots dot --N50 --tsv_stats
  else
    NanoPlot -t 8 --fastq ${sample} --prefix ${sampleid}_raw_ --plots dot --N50 --tsv_stats
  fi
  """
}

process BLASTN {
  tag "${sampleid}"
  containerOptions "${bindOptions}"
  label "setting_10"

  input:
    tuple val(sampleid), path(assembly)
  output:
    tuple val(sampleid), path("${sampleid}*_blastn.bls"), emit: blast_results

  script:
  def blastoutput = assembly.getBaseName() + "_blastn.bls"
  
  if (params.blast_mode == "ncbi") {
    """
    cp ${blastn_db_dir}/taxdb.btd .
    cp ${blastn_db_dir}/taxdb.bti .
    blastn -query ${assembly} \
      -db ${params.blastn_db} \
      -out ${blastoutput} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames' \
      -max_target_seqs 1
    """
  }
  
  else if (params.blast_mode == "localdb") {
    """
    blastn -query ${assembly} \
      -db ${params.blastn_db} \
      -out ${blastoutput} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
      -max_target_seqs 1
    """
  }
}
