process NANOPLOT {
  publishDir { "${params.outdir}/${sampleid}/qc/nanoplot" },  pattern: '{*NanoPlot-report.html}', mode: 'link'
  publishDir { "${params.outdir}/${sampleid}/qc/nanoplot" },  pattern: '{*NanoStats.txt}', mode: 'link'
  publishDir { "${params.outdir}/${sampleid}/qc/nanoplot" },  pattern: '{*LengthvsQualityScatterPlot_dot.html}', mode: 'link'
  tag "${sampleid}"
  label "setting_10"

  input:
    tuple val(sampleid), path(sample)
  output:
    path("*NanoPlot-report.html")
    path("*NanoStats.txt")
    path("*LengthvsQualityScatterPlot_dot.html")
    path("*NanoStats.txt"), emit: read_counts
    tuple val(sampleid), path("${sampleid}_raw_NanoStats.txt"), emit: stats, optional: true

  
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