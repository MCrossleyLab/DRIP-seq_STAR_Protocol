#!/usr/bin/env nextflow
nextflow.enable.dsl=2

process MACS2_CALLPEAK {
    tag "${contrast_name}"
    
    publishDir "${params.outdir}/callpeaks_macs2/${contrast_name}", mode: 'copy'

    input:
    tuple val(contrast_name), path(ip_bams), path(input_bams)

    output:
    path "${contrast_name}_peaks.broadPeak", emit: broadpeaks
    path "${contrast_name}_broadpeaks.log", emit: log
    path "${contrast_name}*", emit: all_outputs

    script:
    """
    samtools merge _merged_ip.bam ${ip_bams.join(' ')}
    samtools merge _merged_in.bam ${input_bams.join(' ')}
    
    macs2 callpeak \\
        --broad \\
        -g mm \\
        -f BAMPE \\
        -t _merged_ip.bam \\
        -c _merged_in.bam \\
        --outdir . \\
        -n ${contrast_name} \\
        2> ${contrast_name}_broadpeaks.log
    """
}
