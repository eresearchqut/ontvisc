if (params.blastn_db != null) {
    blastn_db_name = file(params.blastn_db).name
    blastn_db_dir = file(params.blastn_db).parent
}

switch (workflow.containerEngine) {
  case "singularity":
    bindbuild = "";
    if (params.blastn_db != null) {
      bindbuild = (bindbuild + "-B ${blastn_db_dir} ")
    }
    bindOptions = bindbuild;
    break;
  default:
    bindOptions = "";
}

process MINIMAP2_ALIGN {
  tag "${sampleid}"
  label "setting_8"

  input:
  tuple val(sampleid), path(fastq)
  path(reference)
  output:
  tuple val(sampleid), path(fastq), path("${sampleid}_unaligned_ids.txt"), emit: sequencing_ids
  //tuple val(sampleid), path("${sampleid}.sam"), emit: sam

  def minimap_options = (params.minimap_options) ? " ${params.minimap_options}" : ''

  script:
  """
  minimap2 -ax ${minimap_options} -uf -k14 ${reference} ${fastq} -t ${task.cpus} > ${sampleid}.sam
  awk '\$6 == "*" { print \$0 }' ${sampleid}.sam | cut -f1 | uniq >  ${sampleid}_unaligned_ids.txt
  """
}
/*
tuple val(sampleid), path(fastq), path("${sampleid}_unaligned_ids.txt"), emit: sequencing_ids
awk '\$6 == "*" { print \$0 }' ${sampleid}.sam | cut -f1 | uniq >  ${sampleid}_unaligned_ids.txt
*/


process EXTRACT_READS {
  tag "${sampleid}"
  label "setting_2"

  input:
  tuple val(sampleid), path(fastq), path(unaligned_ids)
  output:
  path("*reads_count.txt"), emit: read_counts
  file("*reads_count.txt")
  tuple val(sampleid), path("*_unaligned.fastq"), emit: unaligned_fq

  script:
  """
  seqtk subseq ${fastq} ${unaligned_ids} > ${sampleid}_unaligned.fastq
  
  n_lines=\$(expr \$(cat ${sampleid}_unaligned.fastq | wc -l) / 4)
  echo \$n_lines > ${sampleid}_unaligned_reads_count.txt
  """
}

/*
n_lines=\$(bc <<< \$(cat ${fastq.baseName}_unaligned.fastq | wc -l)/4)
process MAP_BACK_TO_ASSEMBLY {
  cpus "${params.minimap2_threads}"
  tag "${sampleid}"
  label "large"
  publishDir "$params.outdir/$sampleid",  mode: 'copy', pattern: '*.sam'

  container 'quay.io/biocontainers/minimap2:2.24--h7132678_1'

  input:
  tuple val(sampleid), path(fastq), path(assembly)
  output:
  path "${sampleid}.sam"
  tuple val(sampleid), path(fastq), path("${sampleid}_unaligned_ids.txt"), emit: unmapped_ids

  script:
  """
  minimap2 -ax splice -uf -k14 ${assembly} ${fastq} > ${sampleid}.sam
  awk '\$6 == "*" { print \$0 }' ${sampleid}.sam | cut -f1 | uniq >  ${sampleid}_unaligned_ids.txt
  """
}
*/

process FASTQ2FASTA {
  tag "${sampleid}"
  label "setting_2"
  //publishDir "$params.outdir/$sampleid/fasta", mode: 'copy', pattern: '*.fasta'

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
  publishDir "${params.outdir}/${sampleid}/nanoplot",  pattern: '{*.html,*.txt}', mode: 'link', saveAs: { filename -> "${sampleid}_$filename" }
//  publishDir "${params.outdir}/${sampleid}/nanoplot",  pattern: '*.NanoStats.txt', mode: 'link', saveAs: { filename -> "${sampleid}_$filename" }
  tag "${sampleid}"
  label "setting_2"

  input:
    tuple val(sampleid), path(sample)
  output:
    path("*.html")
    path("*NanoStats.txt")
    path("*NanoStats.txt"), emit: read_counts
  
  script:
 // if (sample.endsWith("quality_trimmed.fastq.gz")) {
  """
  if [[ ${sample} == *trimmed.fastq.gz ]];
  then
    NanoPlot -t 2 --fastq ${sample} --prefix ${sampleid}_filtered_ --plots dot --N50 --tsv_stats
    
  else
    NanoPlot -t 2 --fastq ${sample} --prefix ${sampleid}_raw_ --plots dot --N50 --tsv_stats
  fi
  
  """
}
/*
\$(expr \$(zcat ${sample}  | wc -l) / 4) >> ${sampleid}_quality_trimmed_reads_count.txt
\$(expr \$(zcat ${sample}  | wc -l) / 4) >> ${sampleid}_raw_reads_count.txt
  fi
*/
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