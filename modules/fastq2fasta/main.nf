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