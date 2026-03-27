#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * 1) Count non-duplicate mapped reads per strand
 */
process COUNT_READS_RPM {
    tag "${sample_id}-${species}-${treatment}-${replicate}"

    input:
    tuple val(sample_id), val(species), path(bam_file), val(treatment), val(replicate)

    output:
    tuple val(sample_id), val(species), val(treatment), val(replicate), env(non_duplicate_reads)

    script:
    """
    set -euo pipefail

    # non-duplicate mapped reads:
    #  -F 4   : exclude unmapped
    #  -F 1024: exclude duplicates
    non_duplicate_reads=\$(samtools view -c -F 4 -F 1024 "${bam_file}")

    echo "\$non_duplicate_reads"
    """
}

/*
 * 4) Make bigWig directly from BAM using deepTools bamCoverage
 *    One bigWig per sample / species / strand / treatment / replicate.
 */
process BAMCOVERAGE {
    tag "${sample_id}-${species}-${treatment}-${replicate}"

    publishDir "${params.outdir}/normalized_BigWig_files", mode: 'copy', pattern: "*.bamCov.CPM.bw"

    input:
    tuple val(sample_id), val(species), path(bam_file), val(treatment), val(replicate), val(non_duplicate_reads), val(rpm_scale_factor)
    // tuple val(sample_id), val(species), path(bam_file), val(treatment), val(replicate)

    output:
    tuple val(sample_id), val(species), val(treatment), val(replicate),
          path("${sample_id}.${species}.${treatment}.${replicate}.bamCov.CPM.bw")

    script:
    """
    set -euo pipefail

    # Make sure BAM is indexed (bamCoverage needs an index)
    if [ ! -f "${bam_file}.bai" ]; then
      samtools index "${bam_file}"
    fi

    bamCoverage \\
      --bam "${bam_file}" \\
      --outFileName "${sample_id}.${species}.${treatment}.${replicate}.bamCov.CPM.bw" \\
      --numberOfProcessors max \\
      --normalizeUsing CPM \\
      --outFileFormat bigwig \\
      --binSize 1
    """
}
