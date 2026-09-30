process BLASTN {
  tag "${sampleid}"
  label "setting_10"
  containerOptions "--bind ${file(params.blastn_db).parent}"

  input:
    tuple val(sampleid), path(assembly)
    val(db)

  output:
    tuple val(sampleid), path("${sampleid}*_blastn.bls"), emit: blast_results

  script:
  def blastoutput = assembly.getBaseName() + "_blastn.bls"
  def blastdb_dir  = file(db).parent
  def blastdb_name = file(db).name
  
  if (params.blast_mode == "ncbi") {
    """
    export BLASTDB=${blastdb_dir}
    #cp ${blastdb_dir}/taxdb.btd .
    #cp ${blastdb_dir}/taxdb.bti .
    blastn -query ${assembly} \
      -db ${blastdb_name} \
      -out ${blastoutput} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe sscinames' \
      -max_target_seqs 1
    """
  }
  
  else if (params.blast_mode == "localdb") {
    """
    export BLASTDB=${blastdb_dir}
    blastn -query ${assembly} \
      -db ${blastdb_name}  \
      -out ${blastoutput} \
      -evalue 1e-3 \
      -num_threads ${params.blast_threads} \
      -outfmt '6 qseqid sgi sacc length pident mismatch gapopen qstart qend qlen sstart send slen sstrand evalue bitscore qcovhsp stitle staxids qseq sseq sseqid qcovs qframe sframe' \
      -max_target_seqs 1
    """
  }
}
