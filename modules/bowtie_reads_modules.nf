#!/usr/bin/env nextflow
nextflow.enable.dsl=2

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
