#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

def helpMessage () {
    log.info """
    ONTViSc - ONT-based Viral Screening for Biosecurity
    Marie-Emilie Gauthier 23/05/2023
    Magda Antczak
    Craig Windell
    Roberto Barrero

    Usage:
    Run the command
    nextflow run eresearchqut/ontvisc {optional arguments}...

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

      #### Pre-processing and QC options ####
      --qc_only                       Only perform preliminary QC step using Nanoplot
      --analysis_mode                 clustering, denovo_assembly, read_classification, map2ref
                                      Default: ''
      --adapter_trimming              Run porechop step
                                      [False]
      --porechop_options              Porechop_ABI options
                                      [null]
      --qual_filt                     Run quality filtering step
                                      [False]
      --chopper_options               Chopper options
                                      [null]
      --host_filtering                Run host filtering step using Minimap2
                                      Default: false
      --host_fasta                    Fasta file of nucleotide sequences to filter
                                      [null]
      --canu                          Use Canu for de novo assembly step
                                      Default: false
      --canu_options                  Canu options
                                      Default:  ''
      --flye                          Use Flye for de novo assembly step
      --flye_ont_mode                 Select from nano-raw, nano-corr, nano-hq
                                      Default:  'nano-raw'
      --flye_options                  Flye options
                                      Default: ''
      --rattle_clustering_options     Rattle clustering options
                                      Default: ''
      --rattle_polishing_options      Rattle polishing options
                                      Default: ''
      --final_primer_check            Performs a final primer check after performing de novo assembly step
                                      Default: false

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
if (params.host_fasta != null) {
    host_fasta_dir = file(params.host_fasta).parent
}
if (params.kaiju_nodes != null & params.kaiju_dbname != null & params.kaiju_names != null) {
    kaiju_dbs_dir = file(params.kaiju_nodes).parent
}
if (params.krkdb != null) {
    krkdb_dir = file(params.krkdb).parent
}

if (params.porechop_custom_primers == true) {
    porechop_custom_primers_dir = file(params.porechop_custom_primers_path).parent
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
    if (params.host_fasta != null) {
      bindbuild = (bindbuild + "-B ${host_fasta_dir} ")
    }
    if (params.porechop_custom_primers) {
      bindbuild = (bindbuild + "-B ${porechop_custom_primers_dir} ")
    }
    bindOptions = bindbuild;
    break;
  default:
    bindOptions = "";
}
/*
process MERGE {
  //publishDir "${params.outdir}/${sampleid}/merge", pattern: '*.fastq.gz', mode: 'link'
  tag "${sampleid}"
  label 'setting_1'

  input:
    tuple val(sampleid), path(lanes)
  output:
    tuple val(sampleid), path("${sampleid}.fastq.gz"), emit: merged
  script:
  """
  cat ${lanes} > ${sampleid}.fastq.gz

  """
}
*/
process QCREPORT {
  publishDir "${params.outdir}/qc_report", mode: 'copy', overwrite: true
  containerOptions "${bindOptions}"

  input:
    path multiqc_files

  output:
    path("run_qc_report_*txt")
    path("run_qc_report_*html")

  script:
    """
    seq_run_qc_report.py --host_filtering ${params.host_filtering} --adapter_trimming ${params.adapter_trimming} --quality_trimming ${params.qual_filt}
    """
}

process CUTADAPT {
  publishDir "${params.outdir}/${sampleid}/assembly", pattern: '*_filtered.fa', mode: 'link'
  publishDir "${params.outdir}/${sampleid}/assembly", pattern: '*_cutadapt.log', mode: 'link'
  tag "${sampleid}"
  label 'medium'

  input:
    tuple val(sampleid), path(sample)

  output:
    file("${sampleid}_cutadapt.log")
    file("${sample.baseName}_filtered.fa")
    tuple val(sampleid), path("${sample.baseName}_filtered.fa"), emit: trimmed

  script:
    """
    cutadapt -j ${task.cpus} -g "AAGCAGTGGTATCAACGCAGAGTACGCGGG;min_overlap=14" -a "CCCGCGTACTCTGCGTTGATACCACTGCTT;min_overlap=14" -o ${sample.baseName}_filtered.fa ${sample} > ${sampleid}_cutadapt.log
    sed -i 's/len=[0-9]* reads/reads/' ${sample.baseName}_filtered.fa
    sed -i 's/ trim=.*\$//' ${sample.baseName}_filtered.fa
    """
}

process CHOPPER {
  publishDir "${params.outdir}/${sampleid}/preprocessing/chopper", pattern: '*_chopper.log', mode: 'link'
  tag "${sampleid}"
  label 'setting_3'

  input:
    tuple val(sampleid), path(sample)

  output:
    path("${sampleid}_chopper.log")
    tuple val(sampleid), path("${sampleid}_filtered.fastq.gz"), emit: chopper_filtered_fq

  script:
  def chopper_options = (params.chopper_options) ? " ${params.chopper_options}" : ''
    """
    gunzip -c ${sample} | chopper ${chopper_options} 2> ${sampleid}_chopper.log | gzip > ${sampleid}_filtered.fastq.gz
    """
}

/*
#recommended settings for CANU using metagnomics data
https://github.com/marbl/canu/issues/2079

genomeSize=20m
maxInputCoverage=10000 corOutCoverage=10000
corMhapSensitivity=high
corMinCoverage=0
redMemory=32 oeaMemory=32 batMemory=64
useGrid=false
minReadLength=200
minOverlapLength=50
maxThreads=4
minInputCoverage=0
stopOnLowCoverage=0

maxInputCoverage=10000 corOutCoverage=10000 corMhapSensitivity=high corMinCoverage=0 redMemory=32 oeaMemory=32 batMemory=64 useGrid=false minReadLength=200 minOverlapLength=50 maxThreads=4 minInputCoverage=0 stopOnLowCoverage=0
*/
process CANU {
  publishDir "${params.outdir}/${sampleid}/assembly/canu", mode: 'copy', pattern: '{*.fasta,*.log}'
  tag "${sampleid}"
  label 'setting_8'

  input:
    tuple val(sampleid), path(fastq)

  output:
    path("${sampleid}_canu_assembly.fasta")
    path("${sampleid}_canu.log")
    tuple val(sampleid), path("${sampleid}_canu.fastq.gz"), path("${sampleid}_canu_assembly.fasta"), emit: assembly
    tuple val(sampleid), path("${sampleid}_canu_assembly.fasta"), emit: assembly2


  script:
  def canu_options = (params.canu_options) ? " ${params.canu_options}" : ''
    """
    canu -p ${sampleid} -d ${sampleid} \
      genomeSize=${params.canu_genome_size} \
      -nanopore ${fastq} ${canu_options} 2> ${sampleid}_canu.log

    if [[ ! -s ${sampleid}/${sampleid}.contigs.fasta ]]
      then
        touch ${sampleid}_canu_assembly.fasta
    else
      cat ${sampleid}/${sampleid}.contigs.fasta ${sampleid}/${sampleid}.unassembled.fasta > ${sampleid}_canu_assembly.fasta
    fi
    cp ${fastq} ${sampleid}_canu.fastq.gz
    """
}

process FLYE {
  publishDir "${params.outdir}/${sampleid}/assembly/flye", mode: 'copy', pattern: '{*.fasta,*.log}'
  tag "${sampleid}"
  label 'setting_9'

  input:
    tuple val(sampleid), path(fastq)

  output:
    path("${sampleid}_flye_assembly.fasta")
    path("${sampleid}_flye.log")
    tuple val(sampleid), path("${sampleid}_flye.fastq.gz"), path("${sampleid}_flye_assembly.fasta"), emit: assembly
    tuple val(sampleid), path("${sampleid}_flye_assembly.fasta"), emit: assembly2

  script:
  def flye_options = (params.flye_options) ? " ${params.flye_options}" : ''
    """
    flye ${flye_options} --out-dir outdir --threads ${task.cpus} --${params.flye_mode} ${fastq}

    if [[ ! -s outdir/assembly.fasta ]]
      then
        touch ${sampleid}_flye_assembly.fasta
    else
      cp outdir/assembly.fasta ${sampleid}_flye_assembly.fasta
    cp outdir/flye.log ${sampleid}_flye.log
    fi
    cp ${fastq} ${sampleid}_flye.fastq.gz
    """
}
/*
errorStrategy 'ignore'
*/

process BLASTN2REF {
  publishDir "${params.outdir}/${sampleid}", mode: 'copy', pattern: '*/*/*txt'
  tag "${sampleid}"
  label 'setting_1'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(assembly)

  output:
    path "*/*/${sampleid}_blastn_reference_vs_*.txt"

  script:
    """
    if [[ ${assembly} == *_assembly*.fa* ]] ;
    then
      if [[ ${assembly} == *canu_assembly*.fa* ]] ;
      then
        mkdir -p assembly/blast_to_ref
        blastn -query ${assembly} -subject ${reference_dir}/${reference_name} -evalue 1e-3 -out ${sampleid}_blastn_reference_vs_canu_assembly_tmp.txt \
        -outfmt '6 qseqid sacc length pident mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp qcovs' -max_target_seqs 5

        echo "qseqid\tsacc\tlength\tpident\tmismatch\tgapopen\tqstart\tqend\tqlen\tsstart\tsend\tslen\tevalue\tbitscore\tqcovhsp\tqcovs" > header

        cat header ${sampleid}_blastn_reference_vs_canu_assembly_tmp.txt >  assembly/blast_to_ref/${sampleid}_blastn_reference_vs_canu_assembly.txt
      elif [[ ${assembly} == *flye_assembly*.fasta ]] ;
      then
        mkdir -p assembly/blast_to_ref

        blastn -query ${assembly} -subject ${reference_dir}/${reference_name} -evalue 1e-3 -out ${sampleid}_blastn_reference_vs_flye_assembly_tmp.txt \
        -outfmt '6 qseqid sacc length pident mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp qcovs' -max_target_seqs 5

        echo "qseqid\tsacc\tlength\tpident\tmismatch\tgapopen\tqstart\tqend\tqlen\tsstart\tsend\tslen\tevalue\tbitscore\tqcovhsp\tqcovs" > header

        cat header ${sampleid}_blastn_reference_vs_flye_assembly_tmp.txt > assembly/blast_to_ref/${sampleid}_blastn_reference_vs_flye_assembly.txt
      fi
    elif [[ ${assembly} == *clustering.fasta ]] ;
    then
      mkdir -p clustering/blast_to_ref

      blastn -query ${assembly} -subject ${reference_dir}/${reference_name} -evalue 1e-3 -out ${sampleid}_blastn_reference_vs_clustering_tmp.txt \
      -outfmt '6 qseqid sacc length pident mismatch gapopen qstart qend qlen sstart send slen evalue bitscore qcovhsp qcovs' -max_target_seqs 5

      echo "qseqid\tsacc\tlength\tpident\tmismatch\tgapopen\tqstart\tqend\tqlen\tsstart\tsend\tslen\tevalue\tbitscore\tqcovhsp\tqcovs" > header

      cat header ${sampleid}_blastn_reference_vs_clustering_tmp.txt > clustering/blast_to_ref/${sampleid}_blastn_reference_vs_clustering_assembly.txt
    fi
    """
}

process MINIMAP2_REF {
  tag "${sampleid}"
  label 'setting_2'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(sample)

  output:
    tuple val(sampleid), file("${sampleid}_aln.sam"), emit: aligned_sample

  script:
    """
    minimap2 -ax map-ont --MD --sam-hit-only ${reference_dir}/${reference_name} ${sample} > ${sampleid}_aln.sam
    """
}

process SAMTOOLS {
  publishDir "${params.outdir}/${sampleid}/mapping", mode: 'copy'
  tag "${sampleid}"
  label 'setting_2'

  input:
    tuple val(sampleid), path(sample)

  output:
    path "${sampleid}_aln.sorted.bam"
    path "${sampleid}_aln.sorted.bam.bai"
    path "${sampleid}_coverage.txt"
    path "${sampleid}_histogram"
    tuple val(sampleid), path("${sampleid}_aln.sorted.bam"), path("${sampleid}_aln.sorted.bam.bai"), emit: sorted_sample

  script:
    """
    samtools view -Sb -F 4 ${sample} | samtools sort -o ${sampleid}_aln.sorted.bam
    samtools index ${sampleid}_aln.sorted.bam
    samtools coverage ${sampleid}_aln.sorted.bam > ${sampleid}_histogram.txt  > ${sampleid}_coverage.txt
    samtools coverage -A -w 50 ${sampleid}_aln.sorted.bam > ${sampleid}_histogram
    """
}

process EXTRACT_REF_FASTA {
  tag "$sampleid"
  label "setting_1"
  publishDir "${params.outdir}/${sampleid}/alignments", mode: 'copy', pattern: '*fasta'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results)

  output:
    path("*fasta"), optional: true
    tuple val(sampleid), path("*fasta"), emit: fasta_files, optional: true
  
  script:
    """
    cut -f1,4 ${blast_results} | sed '1d' | sed 's/ /_/g' > ids_to_retrieve.txt
    if [ -s ids_to_retrieve.txt ]
      then
        for i in `cut -f2  ids_to_retrieve.txt`; do j=`grep \${i} ids_to_retrieve.txt | cut -f1`; efetch -db nucleotide  -id \${i} -format fasta > ${sampleid}_\${i}_\${j}.fasta ; done
    fi
    """
}

process MAPPING_BACK_TO_REF {
  tag "$sampleid"
  label "setting_3"
  publishDir "${params.outdir}/${sampleid}/alignments", mode: 'copy', pattern: '*sorted.bam*'
  //publishDir "${params.outdir}/01_VirReport/${sampleid}/alignments/NT", mode: 'link', overwrite: true, pattern: "*{.fa*,.fasta,metrics.txt,scores.txt,targets.txt,stats.txt,log.txt,.bcf*,.vcf.gz*,.bam*}"

  input:
    tuple val(sampleid), path(results)

  output:
    path("*bam"), optional: true
    path("*bam.bai"), optional: true
    tuple val(sampleid), path("*sorted.bam"), emit: bam_files, optional: true
    tuple val(sampleid), path("*sorted.bam.bai"), emit: bai_files, optional: true

  script:
    """
    if compgen -G "*.fasta" > /dev/null;
      then
        mapping_back_to_ref.py --fastq ${sampleid}_preprocessed.fastq.gz
    fi
    """
}

process MOSDEPTH {
  tag "$sampleid"
  label "setting_3"
  publishDir "${params.outdir}/${sampleid}/alignments", mode: 'copy'

  input:
    tuple val(sampleid), path("*")

  output:
    path("*"), optional: true
    tuple val(sampleid), path("*mosdepth.global.dist.txt"), emit: mosdepth_results, optional: true

  script:
    """
    if compgen -G "*.bam" > /dev/null;
      then
        for i in *bam;
        do
          echo \${i%.sorted.bam}
          filen=`echo "\${i%.sorted.bam}"`
          mosdepth \${filen} \${i};
        done
    fi
    """
}

process COVERM {
  tag "$sampleid"
  label "setting_3"
  publishDir "${params.outdir}/${sampleid}/alignments", mode: 'copy'

  input:
    tuple val(sampleid), path("*")

  output:
    path("*"), optional: true
    tuple val(sampleid), path("*_coverm_summary.txt"), emit: coverm_results, optional: true

  script:
    """
    if compgen -G "*.bam" > /dev/null;
      then
        for i in *bam;
        do
          echo \${i%.sorted.bam}
          filen=`echo "\${i%.sorted.bam}"`
          coverm genome --genome-fasta-files \${filen}.fasta --bam-files \${i} --threads ${task.cpus} --output-file \${filen}_coverm_summary.txt -m count mean variance rpkm covered_bases length --min-covered-fraction 0;
          coverm genome --genome-fasta-files \${filen}.fasta --bam-files \${i} --threads ${task.cpus} --output-file \${filen}_coverage_histogram.txt -m coverage_histogram --min-covered-fraction 0;
        done
    fi
    """
}

process COVSTATS {
  tag "$sampleid"
  label "setting_2"
  publishDir "${params.outdir}/${sampleid}/alignments", mode: 'copy'

  input:
    tuple val(sampleid), path("*")

  output:
    path("*top_blast_with_cov_stats.txt"), optional: true
    path("*top_blast_with_cov_stats.txt"), emit: detections_summary, optional: true

  script:
    """
    if compgen -G "*coverm_summary.txt" > /dev/null;
      then
        derive_coverage_stats.py --sample ${sampleid}
    fi
    """
}

/*
process BAMCOVERAGE {
  publishDir "${params.outdir}/${sampleid}/mapping", mode: 'link'
  tag "${sampleid}"
  label 'setting_11'

  input:
    tuple val(sampleid), path(bam), path(bai)
  output:
    path "${sampleid}.bw"
  script:
  """
  bamCoverage -b ${bam} -o ${sampleid}.bw
  """
}


process BAMCOVERAGE {
  publishDir "${params.outdir}/${sampleid}/mapping", mode: 'link'
  tag "${sampleid}"
  label 'setting_11'

  input:
    tuple val(sampleid), path(bam), path(bai)
  output:
    path "${sampleid}.bamcov.txt"
  script:
  """
  bamcov -H ${bam} > ${sampleid}.bamcov.txt
  """
}
*/

process PORECHOP_ABI {
  tag "${sampleid}"
  publishDir "$params.outdir/${sampleid}/preprocessing/porechop",  mode: 'link', pattern: '*_porechop.log'
  label "setting_9"
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(sample)

  output:
    file("${sampleid}_porechop_trimmed.fastq.gz")
    file("${sampleid}_porechop.log")
    tuple val(sampleid), file("${sampleid}_porechop_trimmed.fastq.gz"), emit: porechopabi_trimmed_fq

  script:
  def porechop_options = (params.porechop_options) ? " ${params.porechop_options}" : ''
    """
    if [[ ${params.porechop_custom_primers} == true ]]; then
      porechop_abi -i ${sample} -t ${task.cpus} -o ${sampleid}_porechop_trimmed.fastq.gz --custom_adapters ${params.porechop_custom_primers_path} ${porechop_options}  > ${sampleid}_porechop.log
    else
      porechop_abi -i ${sample} -t ${task.cpus} -o ${sampleid}_porechop_trimmed.fastq.gz ${porechop_options}  > ${sampleid}_porechop.log
    fi
    """
}

//trims fastq read names after the first whitespace
process REFORMAT {
  tag "${sampleid}"
  label "setting_3"
  publishDir "$params.outdir/${sampleid}/preprocessing", mode: 'copy'

  input:
    tuple val(sampleid), path(fastq)

  output:
    tuple val(sampleid), path("${sampleid}_preprocessed.fastq.gz"), emit: reformatted_fq
    tuple val(sampleid), path("${sampleid}_preprocessed.fastq.gz"), emit: cov_derivation_ch

  script:
    """
    reformat.sh in=${fastq} out=${sampleid}_preprocessed.fastq.gz trd qin=33
    """
}

process CAP3 {
  tag "${sampleid}"
  label "setting_3"
  time "3h"
  publishDir "$params.outdir/$sampleid/clustering/cap3", mode: 'copy', pattern: '*_clustering.fasta'
  publishDir "$params.outdir/$sampleid/clustering/rattle", mode: 'copy', pattern: '*_rattle.fasta'

  input:
    tuple val(sampleid), path(fasta)

  output:
    file("${sampleid}_rattle.fasta")
    tuple val(sampleid), path("${sampleid}_clustering.fasta"), emit: scaffolds

  script:
    """
    cap3 ${fasta}
    cat ${fasta}.cap.singlets ${fasta}.cap.contigs > ${sampleid}_clustering.fasta
    cp ${fasta} ${sampleid}_rattle.fasta
    """
}

process EXTRACT_VIRAL_BLAST_HITS {
  tag "${sampleid}"
  label "setting_2"
  publishDir "$params.outdir/$sampleid", overwrite: true
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(blast_results)

  output:
    file "*/*/${sampleid}*_blastn_top_hits.txt"
    file "*/*/${sampleid}*_blastn_top_viral_hits*.txt"
    file "*/*/${sampleid}*_blastn_top_viral_spp_hits.txt"
    file "*/*/${sampleid}*_queryid_list_with_viral_match.txt"
    file "*/*/${sampleid}*_viral_spp_abundance*.txt"
    file "*/*/*report*html"

    tuple val(sampleid), path("*/*/${sampleid}*_blastn_top_viral_spp_hits.txt"), path("*/*/${sampleid}*_queryid_list_with_viral_match.txt"), path("*/*/${sampleid}*_viral_spp_abundance.txt"), emit: blast_results
    tuple val(sampleid), path("*/*/${sampleid}*_blastn_top_viral_spp_hits.txt"), emit: blast_results2

  script:
    """
    if [[ ${params.analysis_mode} == "clustering" ]] ;
    then
      mkdir -p clustering/blastn
      cat ${blast_results} > clustering/blastn/${sampleid}_blastn.txt
      cd clustering/blastn
      select_top_blast_hit.py --sample_name ${sampleid} --blastn_results ${sampleid}_blastn.txt --analysis_method clustering --mode ${params.blast_mode}
    elif [[ ${params.analysis_mode} == "denovo_assembly" ]] ;
    then
      mkdir -p assembly/blastn
      cat ${blast_results} > assembly/blastn/${sampleid}_blastn.txt
      cd assembly/blastn
      select_top_blast_hit.py --sample_name ${sampleid} --blastn_results ${sampleid}_blastn.txt --analysis_method assembly --mode ${params.blast_mode}
    elif [[ ${params.analysis_mode} == "read_classification" ]]
      mkdir -p read_classification/homology_search
      cat ${blast_results} > read_classification/homology_search/${sampleid}_blastn.txt
      cd read_classification/homology_search
      select_top_blast_hit.py --sample_name ${sampleid} --blastn_results ${sampleid}_blastn.txt --analysis_method read_classification --mode ${params.blast_mode}
    fi
    """
}

process CONCATENATE_FASTA {
  tag "${sampleid}"
  label "setting_2"
  publishDir "${params.outdir}/${sampleid}", mode: 'link'

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
  publishDir "${params.outdir}/${sampleid}/read_classification/kaiju", mode: 'copy'
  label 'setting_4'
  containerOptions "${bindOptions}"
  tag "${sampleid}"

  input:
    tuple val(sampleid), path(fastq)

  output:
    file "${sampleid}_kaiju_name.tsv"
    file "${sampleid}_kaiju_summary*.tsv"
    file "${sampleid}_kaiju.krona"
    tuple val(sampleid), path("${sampleid}_kaiju_summary_viral.tsv"), emit: kaiju_results
    tuple val(sampleid), path("*kaiju.krona"), emit: krona_results

  script:
    """
    c1grep() { grep "\$@" || test \$? = 1; }

    kaiju \
        -z ${params.kaiju_threads} \
        -t ${params.kaiju_nodes}  \
        -f ${params.kaiju_dbname} \
        -o ${sampleid}_kaiju.tsv \
        -i ${fastq} \
        -v

    kaiju-addTaxonNames -t ${params.kaiju_nodes} -n ${params.kaiju_names} -i ${sampleid}_kaiju.tsv -o ${sampleid}_kaiju_name.tsv
    kaiju2table -e -t ${params.kaiju_nodes} -r species -n ${params.kaiju_names} -o ${sampleid}_kaiju_summary.tsv ${sampleid}_kaiju.tsv
    kaiju2krona -t ${params.kaiju_nodes} -n ${params.kaiju_names} -i ${sampleid}_kaiju.tsv -o ${sampleid}_kaiju.krona

    c1grep "taxon_id\\|virus\\|viroid\\|viricota\\|viridae\\|viriform\\|virales\\|virinae\\|viricetes\\|virae\\|viral" ${sampleid}_kaiju_summary.tsv > ${sampleid}_kaiju_summary_viral.tsv
    awk -F'\\t' '\$2>=0.05' ${sampleid}_kaiju_summary_viral.tsv > ${sampleid}_kaiju_summary_viral_filtered.tsv
    """
}
/*
process KAIJU_HTML {
  tag "${sampleid}"
  label "local"
  publishDir "$params.outdir/$sampleid/read_classification/kaiju",  mode: 'link'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(kaiju_report)

  output:
    file("*_kaiju_report.html")


  script:
  """
  kraken_html_report.py  --sample ${sampleid}
  """
}
*/
process READ_CLASSIFICATION_HTML {
  publishDir "${params.outdir}/${sampleid}/read_classification/summary", mode: 'copy', overwrite: true
  label 'local'
  containerOptions "${bindOptions}"
  tag "${sampleid}"

  input:
    tuple val(sampleid), path(results)

  output:
    file "*_read_classification_report.html"

  script:
    """
    summary_read_classification.py --sample ${sampleid}
    """
}

process KRONA {
  publishDir "${params.outdir}/${sampleid}/read_classification/kaiju", mode: 'link'
  label 'setting_3'
  containerOptions "${bindOptions}"
  tag "${sampleid}"

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
  label 'setting_5'
  publishDir "$params.outdir/$sampleid/read_classification/kraken",  mode: 'copy'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(fastq)

  output:
    file("${sampleid}.kraken2")
    file("${sampleid}_kraken_report.txt")
    file("${sampleid}_unclassified.fastq")
    tuple val(sampleid), path("${sampleid}_kraken_report.txt"), emit: results

  script:
    """
    kraken2 --db ${params.krkdb} \
            --use-names \
            --threads ${task.cpus} \
            --report ${sampleid}_kraken_report.txt \
            --report-minimizer-data \
            --minimum-hit-groups 3 \
            --unclassified-out ${sampleid}_unclassified.fastq \
            ${fastq} > ${sampleid}.kraken2
    """
}

process BRACKEN {
  tag "${sampleid}"
  label 'setting_2'
  publishDir "$params.outdir/$sampleid/read_classification/kraken",  mode: 'copy'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(kraken_report)

  output:
    file("${sampleid}_bracken_report*.txt")

    tuple val(sampleid), path("${sampleid}_bracken_report_viral.txt"), emit: bracken_results

  script:
    """
    c1grep() { grep "\$@" || test \$? = 1; }

    est_abundance.py -i ${kraken_report} \
                    -k ${params.krkdb}/database50mers.kmer_distrib \
                    -t 1 \
                    -l S -o ${sampleid}_bracken_report.txt


    c1grep  "taxonomy_id\\|virus\\|viroid" ${sampleid}_bracken_report.txt > ${sampleid}_bracken_report_viral.txt
    awk -F'\\t'  '\$7>=0.0001'  ${sampleid}_bracken_report_viral.txt > ${sampleid}_bracken_report_viral_filtered.txt
    """
}

process BRACKEN_HTML {
  tag "${sampleid}"
	label "local"
	publishDir "$params.outdir/$sampleid/read_classification/kraken",  mode: 'link'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(braken_report)

  output:
    file("*_bracken_report.html")

  script:
    """
    kraken_html_report.py --sample ${sampleid}
    """
}

process RATTLE {
  publishDir "${params.outdir}/${sampleid}/clustering/rattle", mode: 'copy', pattern: 'transcriptome.fq'
  tag "${sampleid}"
  label 'setting_7'
  containerOptions "${bindOptions}"

  input:
    tuple val(sampleid), path(fastq)

  output:
    file("transcriptome.fq")
    tuple val(sampleid), path("transcriptome.fq"), emit: clusters
    tuple val(sampleid), path(fastq), path("transcriptome.fq"), emit: clusters2

  script:
  def rattle_polishing_options = (params.rattle_polishing_options) ? " ${params.rattle_polishing_options}" : ''
  def rattle_clustering_options = (params.rattle_clustering_options) ? " ${params.rattle_clustering_options}" : ''
    """
    rattle cluster -i ${fastq} -t ${task.cpus} ${rattle_clustering_options}  -o .
    rattle cluster_summary -i ${fastq} -c clusters.out > ${sampleid}_cluster_summary.txt
    mkdir clusters
    rattle extract_clusters -i ${fastq} -c clusters.out -l ${sampleid} -o clusters --fastq
    rattle correct -i ${fastq} -c clusters.out -t ${task.cpus} -l ${sampleid}
    rattle polish -i consensi.fq -t ${task.cpus} --summary ${rattle_polishing_options}
    """
}

process MEDAKA {
  publishDir "${params.outdir}/${sampleid}/mapping", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'
  containerOptions "${bindOptions}"

  input:
   tuple val(sampleid), path(bam), path(bai)

  output:
    tuple val(sampleid), path("${sampleid}_medaka.annotated.unfiltered.vcf"), emit: unfilt_vcf

  script:
  def medaka_consensus_options = (params.medaka_consensus_options) ? " ${params.medaka_consensus_options}" : ''
    """
    medaka consensus ${bam} ${sampleid}_medaka_consensus_probs.hdf \
      ${medaka_consensus_options} --threads ${task.cpus}

    medaka variant ${reference_dir}/${reference_name} ${sampleid}_medaka_consensus_probs.hdf ${sampleid}_medaka.vcf
    medaka tools annotate --dpsp ${sampleid}_medaka.vcf ${reference_dir}/${reference_name} ${bam} \
          ${sampleid}_medaka.annotated.unfiltered.vcf
    """
}

process FILTER_VCF {
  publishDir "${params.outdir}/${sampleid}/mapping", mode: 'copy'
  tag "${sampleid}"
  label 'setting_3'
  containerOptions "${bindOptions}"

  input:
   tuple val(sampleid), path(vcf)

  output:
    path("${sampleid}_medaka.consensus.fasta")
    path("${sampleid}_medaka.annotated.vcf.gz")

  script:
    """
    bcftools reheader ${vcf} -s <(echo '${sampleid}') \
    | bcftools filter \
        -e 'INFO/DP < ${params.bcftools_min_coverage}' \
        -s LOW_DEPTH \
        -Oz -o ${sampleid}_medaka.annotated.vcf.gz

    # create consensus
    bcftools index ${sampleid}_medaka.annotated.vcf.gz
    bcftools consensus -f ${reference_dir}/${reference_name} ${sampleid}_medaka.annotated.vcf.gz \
        -i 'FILTER="PASS"' \
        -o ${sampleid}_medaka.consensus.fasta
    """
}


process FASTCAT {
  publishDir "${params.outdir}/${sampleid}/qc/fastcat", mode: 'copy'
  tag "${sampleid}"
  label "setting_2"

  input:
    tuple val(sampleid), path(fastq)

  output:
    path("${sampleid}_stats.tsv")
    path("histograms/*")
    tuple val(sampleid), path("${sampleid}.fastq.gz"), emit: merged

  script:
    """
    fastcat \
        -s ${sampleid} \
        -f ${sampleid}_stats.tsv \
        --histograms histograms \
        ${fastq} \
        | bgzip > ${sampleid}.fastq.gz
    """
}

process DETECTION_REPORT {
  label "local"
    publishDir "${params.outdir}/summary", mode: 'copy', overwrite: true
    containerOptions "${bindOptions}"

  input:
    path('*')

  output:
    path("summary_detection.txt")

  script:
    """
    detection_summary.py --threshold ${params.contamination_flag}
    """
}

include { MINIMAP2_ALIGN_RNA } from './modules.nf'
include { EXTRACT_READS as EXTRACT_READS_STEP1 } from './modules.nf'
include { EXTRACT_READS as EXTRACT_READS_STEP2 } from './modules.nf'
include { EXTRACT_READS as EXTRACT_READS_STEP3 } from './modules.nf'
include { MINIMAP2_ALIGN_DNA as MAP_BACK_TO_CONTIGS } from './modules.nf'
include { MINIMAP2_ALIGN_DNA as MAP_BACK_TO_CLUSTERS } from './modules.nf'
include { FASTQ2FASTA } from './modules.nf'
include { FASTQ2FASTA as FASTQ2FASTA_STEP1} from './modules.nf'
include { FASTQ2FASTA as FASTQ2FASTA_STEP2} from './modules.nf'
include { NANOPLOT as QC_PRE_DATA_PROCESSING } from './modules.nf'
include { NANOPLOT as QC_POST_DATA_PROCESSING } from './modules.nf'
include { BLASTN as READ_CLASSIFICATION_BLASTN } from './modules.nf'
include { BLASTN as ASSEMBLY_BLASTN } from './modules.nf'

workflow {
  if (params.samplesheet) {
    Channel
      .fromPath(params.samplesheet, checkIfExists: true)
      .splitCsv(header:true)
      .map{ row-> tuple((row.sampleid), file(row.sample_files)) }
      .set{ ch_sample }
  } else { exit 1, "Input samplesheet file not specified!" }

  if ( params.analysis_mode == 'clustering' | params.analysis_mode == 'denovo_assembly' | (params.analysis_mode == 'read_classification' & params.megablast)) {
    if (!params.blast_vs_ref) {
      if ( params.blastn_db == null) {
        error("Please provide the path to a blast database using the parameter --blastn_db.")
      }
    }
    else if (params.blast_vs_ref ) {
      if ( params.reference == null) {
      error("Please provide the path to a reference fasta file with the parameter --reference.")
      }
    }
  }
  else if ( params.analysis_mode == 'map2ref' ) {
    if ( params.reference == null) {
      error("Please provide the path to a reference fasta file with the parameter --reference.")
      }
  }

  if (params.merge) {
    FASTCAT ( ch_sample )
    QC_PRE_DATA_PROCESSING ( FASTCAT.out.merged )
    fq = FASTCAT.out.merged
  }
  else {
    fq = ch_sample
    QC_PRE_DATA_PROCESSING ( fq )
  }

  // Data pre-processing
  // Remove adapters uisng either PORECHOP_ABI or CUTADAPT
  if (!params.qc_only) {
    if (params.adapter_trimming) {
      PORECHOP_ABI ( fq )
      trimmed_fq = PORECHOP_ABI.out.porechopabi_trimmed_fq
    }

    else {
      trimmed_fq = fq
    }

    // Quality filtering of reads using chopper
    if (params.qual_filt) {
      CHOPPER ( trimmed_fq)
      filtered_fq = CHOPPER.out.chopper_filtered_fq
    }
    else { filtered_fq = trimmed_fq
    }

    REFORMAT( filtered_fq )

    if ( params.qual_filt & params.adapter_trimming | !params.qual_filt & params.adapter_trimming | params.qual_filt & !params.adapter_trimming) {
      QC_POST_DATA_PROCESSING ( filtered_fq )
    }

    if (params.host_filtering) {
      if ( params.host_fasta == null) {
        error("Please provide the path to a fasta file of host sequences that need to be filtered with the parameter --host_fasta.")
      }
      else {
        MINIMAP2_ALIGN_RNA ( REFORMAT.out.reformatted_fq, params.host_fasta )
        EXTRACT_READS_STEP1 ( MINIMAP2_ALIGN_RNA.out.sequencing_ids )
        final_fq = EXTRACT_READS_STEP1.out.unaligned_fq
      }
    }
    else {
      final_fq = REFORMAT.out.reformatted_fq
    }

    if ( params.qual_filt & params.host_filtering | params.adapter_trimming & params.host_filtering ) {
      ch_multiqc_files = ch_multiqc_files.mix(QC_PRE_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
      ch_multiqc_files = ch_multiqc_files.mix(EXTRACT_READS_STEP1.out.read_counts.collect().ifEmpty([]))
      ch_multiqc_files = ch_multiqc_files.mix(QC_POST_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
      QCREPORT(ch_multiqc_files.collect())
    }

    else if ( params.host_filtering & !params.adapter_trimming & !params.qual_filt ) {
      ch_multiqc_files = ch_multiqc_files.mix(QC_PRE_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
      ch_multiqc_files = ch_multiqc_files.mix(EXTRACT_READS_STEP1.out.read_counts.collect().ifEmpty([]))
      QCREPORT(ch_multiqc_files.collect())
    }

    else if ( params.qual_filt & !params.host_filtering | params.adapter_trimming & !params.host_filtering) {
      ch_multiqc_files = ch_multiqc_files.mix(QC_PRE_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
      ch_multiqc_files = ch_multiqc_files.mix(QC_POST_DATA_PROCESSING.out.read_counts.collect().ifEmpty([]))
      QCREPORT(ch_multiqc_files.collect())
    }

    if (!params.preprocessing_only) {
      //perform clustering using Rattle
      if ( params.analysis_mode == 'clustering' ) {
        RATTLE ( final_fq )
        FASTQ2FASTA( RATTLE.out.clusters )
        CAP3( FASTQ2FASTA.out.fasta )
        contigs = CAP3.out.scaffolds

        if (params.blast_vs_ref) {
          BLASTN2REF ( contigs )
          }
        /*
        else {
          CLUSTERING_BLASTN ( contigs )
          EXTRACT_VIRAL_BLAST_HITS ( CLUSTERING_BLASTN.out.blast_results )
        }
        */
      }

      //perform de novo assembly using either canu or flye
      else if ( params.analysis_mode == 'denovo_assembly' ) {
        if (params.canu) {
          CANU ( final_fq )
          contigs = CANU.out.assembly2
        }
        else if (params.flye) {
          FLYE ( final_fq )
          contigs = FLYE.out.assembly2
        }
        if (params.final_primer_check) {
          CUTADAPT ( contigs )
          contigs = CUTADAPT.out.trimmed
        }
        //limit blast homology search to a reference
        if (params.blast_vs_ref) {
          BLASTN2REF ( contigs )
        }
      
      //blast to NCBI
      }
      if ( params.analysis_mode == 'denovo_assembly' | params.analysis_mode == 'clustering' & !params.blast_vs_ref ) {
        ASSEMBLY_BLASTN ( contigs )
        EXTRACT_VIRAL_BLAST_HITS ( ASSEMBLY_BLASTN.out.blast_results )
        EXTRACT_REF_FASTA (EXTRACT_VIRAL_BLAST_HITS.out.blast_results2)

        mapping_ch = EXTRACT_REF_FASTA.out.fasta_files.concat(REFORMAT.out.cov_derivation_ch).groupTuple().map { [it[0], it[1].flatten()] }//.view()
        MAPPING_BACK_TO_REF ( mapping_ch )
        bamf_ch = MAPPING_BACK_TO_REF.out.bam_files.concat(MAPPING_BACK_TO_REF.out.bai_files, EXTRACT_REF_FASTA.out.fasta_files).groupTuple().map { [it[0], it[1].flatten()] }//.view()
        MOSDEPTH (bamf_ch)
        COVERM (bamf_ch)
        cov_stats_summary_ch = MOSDEPTH.out.mosdepth_results.concat(COVERM.out.coverm_results, EXTRACT_REF_FASTA.out.fasta_files, EXTRACT_VIRAL_BLAST_HITS.out.blast_results2).groupTuple().map { [it[0], it[1].flatten()] }//.view()
        COVSTATS(cov_stats_summary_ch)

        DETECTION_REPORT(COVSTATS.out.detections_summary.collect().ifEmpty([]))
      }

      else if ( params.analysis_mode == 'read_classification') {
      //just perform direct read search
        if (params.megablast) {
          FASTQ2FASTA_STEP1( final_fq )
          READ_CLASSIFICATION_BLASTN( FASTQ2FASTA_STEP1.out.fasta.splitFasta(by: 5000, file: true) )
          READ_CLASSIFICATION_BLASTN.out.blast_results
            .groupTuple()
            .set { ch_blastresults }
          EXTRACT_VIRAL_BLAST_HITS( ch_blastresults )
        }
        if (params.kaiju) {
          KAIJU ( final_fq )
          KRONA ( KAIJU.out.krona_results)
        }
        if (params.kraken2) {
          KRAKEN2 ( final_fq )
          BRACKEN ( KRAKEN2.out.results )
        }

        foo_in_ch = Channel.empty()
        if ( params.megablast & !params.kaiju & !params.kraken2 ) {
        READ_CLASSIFICATION_HTML( EXTRACT_VIRAL_BLAST_HITS.out.blast_results ).concat(EXTRACT_VIRAL_BLAST_HITS.out.blast_results_filt).groupTuple().map { [it[0], it[1].flatten()] }
        }
        else if (params.kaiju & !params.megablast & !params.kraken2) {
          READ_CLASSIFICATION_HTML( KAIJU.out.kaiju_results )
        }
        else if (params.kraken2 & !params.megablast & !params.kaiju) {
          READ_CLASSIFICATION_HTML(BRACKEN.out.bracken_results )
        }
        else if (params.megablast & !params.kaiju & params.kraken2) {
          foo_in_ch = EXTRACT_VIRAL_BLAST_HITS.out.blast_results.concat(BRACKEN.out.bracken_results).groupTuple().map { [it[0], it[1].flatten()] }
          READ_CLASSIFICATION_HTML( foo_in_ch )
        }
        else if (params.megablast & params.kaiju & !params.kraken2) {
          foo_in_ch = EXTRACT_VIRAL_BLAST_HITS.out.blast_results.concat(EXTRACT_VIRAL_BLAST_HITS.out.blast_results_filt, KAIJU.out.kaiju_results).groupTuple().map { [it[0], it[1].flatten()] }
          READ_CLASSIFICATION_HTML( foo_in_ch )
        }
        else if (!params.megablast & params.kaiju & params.kraken2) {
          foo_in_ch = KAIJU.out.kaiju_results.concat(EXTRACT_VIRAL_BLAST_HITS.out.blast_results_filt, BRACKEN.out.bracken_results).groupTuple().map { [it[0], it[1].flatten()] }
          READ_CLASSIFICATION_HTML( foo_in_ch )
        }
        else if (params.megablast & params.kaiju & params.kraken2) {
          foo_in_ch = KAIJU.out.kaiju_results.concat(EXTRACT_VIRAL_BLAST_HITS.out.blast_results, EXTRACT_VIRAL_BLAST_HITS.out.blast_results_filt, BRACKEN.out.bracken_results).groupTuple().map { [it[0], it[1].flatten()] }
          READ_CLASSIFICATION_HTML( foo_in_ch )
        }
      }

      //just perform direct alignment
      else if ( params.analysis_mode == 'map2ref') {
        MINIMAP2_REF ( final_fq )
        SAMTOOLS ( MINIMAP2_REF.out.aligned_sample )
        MEDAKA ( SAMTOOLS.out.sorted_sample )
        FILTER_VCF ( MEDAKA.out.unfilt_vcf )
      }
      else {
        error("Analysis mode (read_classification, clustering, denovo_assembly, map2ref) not specified with e.g. '--analysis_mode clustering' or via a detectable config file.")
      }
    }
  }
}
