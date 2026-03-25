#!/usr/bin/env nextflow
nextflow.enable.dsl=2


process FASTQC_RAW {
    tag "${sample_id}"

    publishDir "${params.outdir}/fastqc_raw", mode: 'copy', pattern: "*_fastqc.*"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("${sample_id}_R1_fastqc.html"), path("${sample_id}_R2_fastqc.html")
    path "${sample_id}_R1_fastqc.zip"
    path "${sample_id}_R2_fastqc.zip"

    script:
    """
    fastqc \\
      --threads ${task.cpus ?: 4} \\
      --outdir . \\
      ${read1} ${read2}
      
    # FastQC names outputs from filenames; rename for consistency
    for f in ${read1} ${read2}; do
      base=\$(basename "\$f")
      
      # Strip common FASTQ extensions (.fastq.gz, .fq.gz, .fastq, .fq)
      prefix="\$base"
      prefix="\${prefix%.fastq.gz}"
      prefix="\${prefix%.fq.gz}"
      prefix="\${prefix%.fastq}"
      prefix="\${prefix%.fq}"
      
      # Decide R1 vs R2 based on prefix
      if [[ "\$prefix" == *R1* || "\$prefix" == *_1 || "\$prefix" == *.1 ]]; then
        mv "\${prefix}_fastqc.html" "${sample_id}_R1_fastqc.html"
        mv "\${prefix}_fastqc.zip"  "${sample_id}_R1_fastqc.zip"
      else
        mv "\${prefix}_fastqc.html" "${sample_id}_R2_fastqc.html"
        mv "\${prefix}_fastqc.zip"  "${sample_id}_R2_fastqc.zip"
      fi
    done
    """
}


process CUTADAPT_ADAPTER_TRIM {
    tag "${sample_id}"

    publishDir "${params.outdir}/cutadapt_adapter", mode: 'copy', pattern: "*.fq.gz"
    publishDir "${params.outdir}/cutadapt_adapter_logs", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(sample_id), path(read1), path(read2)

    output:
    tuple val(sample_id), path("adapt_trim.${sample_id}.R1.fq.gz"), path("adapt_trim.${sample_id}.R2.fq.gz"), emit: sid_fastq_files
    path "adapt_trim.${sample_id}.txt"         // summary log

    script:
    """
    cutadapt \\
      -j ${task.cpus ?: 4} \\
      --minimum-length 30 \\
      -a ${params.adapter_fwd} \\
      -A ${params.adapter_rev} \\
      -o adapt_trim.${sample_id}.R1.fq.gz \\
      -p adapt_trim.${sample_id}.R2.fq.gz \\
      ${read1} ${read2} \\
      1> adapt_trim.${sample_id}.txt
    """
}

process CUTADAPT_TAIL_TRIM {
    tag "${sample_id}"

    publishDir "${params.outdir}/cutadapt_tail", mode: 'copy', pattern: "tail_trim.*.fq.gz"
    publishDir "${params.outdir}/cutadapt_tail_logs", mode: 'copy', pattern: "*.txt"

    input:
    tuple val(sample_id), path(in_r1), path(in_r2)

    output:
    tuple val(sample_id), path("tail_trim.${sample_id}.R1.fq.gz"), path("tail_trim.${sample_id}.R2.fq.gz"), emit: sid_fastq_files
    path "tail_trim.${sample_id}.txt"

    script:
    """
    cutadapt \\
      -j ${task.cpus ?: 4} \\
      --minimum-length 30 \\
      -a "C{30}" \\
      -U 15 \\
      -o tail_trim.${sample_id}.R1.fq.gz \\
      -p tail_trim.${sample_id}.R2.fq.gz \\
      ${in_r1} ${in_r2} \\
      1> tail_trim.${sample_id}.txt
    """
}

process BOWTIE2_ALIGN {
    tag "${sample_id}"

    publishDir "${params.outdir}/bowtie2_bam",  mode: 'copy', pattern: "*.bam"
    publishDir "${params.outdir}/bowtie2_logs", mode: 'copy', pattern: "*.log"

    input:
    tuple val(sample_id), path(r1), path(r2), val(bowtie2_params), val(bowtie2_index), val(species_list)

    output:
    tuple val(sample_id), val(species_list), path("${sample_id}.bam"), emit: sid_bam_file
    path "${sample_id}.bowtie2.log"

    script:
    """
    bowtie2 \\
      -p ${task.cpus ?: 16} \\
      ${bowtie2_params} \\
      -x ${bowtie2_index} \\
      -1 ${r1} -2 ${r2} \\
      2> ${sample_id}.bowtie2.log \\
      | samtools view -bS -@ ${task.cpus} -o ${sample_id}.bam -
    """
}
