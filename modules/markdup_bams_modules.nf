#!/usr/bin/env nextflow
nextflow.enable.dsl=2

/*
 * Input:  tuple(sample_id, bam)
 * Output: tuple(sample_id, collate.bam)
 */
process SAMTOOLS_COLLATE {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(input_bam)

    output:
    tuple val(sample_id), path("*.*")

    script:
    """
    samtools collate -@ ${task.cpus ?: 4} \
      ${input_bam} \
      ${sample_id}.collate.bam
    """
}

/*
 * Input:  tuple(sample_id, collate.bam)
 * Output: tuple(sample_id, fixmate.bam)
 */
process SAMTOOLS_FIXMATE {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(collate_bam)

    output:
    tuple val(sample_id), path("*.*")

    script:
    """
    samtools fixmate -@ ${task.cpus ?: 4} -r -m \
      ${collate_bam} \
      ${sample_id}.fixmate.bam
    """
}

/*
 * Input:  tuple(sample_id, fixmate.bam)
 * Output: tuple(sample_id, positionsorted.bam)
 */
process SAMTOOLS_SORT {
    tag "${sample_id}"

    input:
    tuple val(sample_id), path(fixmate_bam)

    output:
    tuple val(sample_id), path("${sample_id}.positionsorted.bam")

    script:
    """
    samtools sort -@ ${task.cpus ?: 4} \
      ${fixmate_bam} \
      -o ${sample_id}.positionsorted.bam
    """
}

/*
 * Input:  tuple(sample_id, positionsorted.bam)
 * Output: tuple(sample_id, markdup.bam)
 */
process SAMTOOLS_MARKDUP_INDEX {
    tag "${sample_id}"

    publishDir "${params.outdir}/markdup_bams", mode: 'copy', pattern: "*.markdup.bam*"

    input:
    tuple val(sample_id), path(pos_sorted_bam)

    output:
    tuple val(sample_id), path("${sample_id}.markdup.bam")

    script:
    """
    samtools markdup -@ ${task.cpus ?: 4} \
      ${pos_sorted_bam} \
      ${sample_id}.markdup.bam
    samtools index ${sample_id}.markdup.bam
    """
}

/*
 * Input:  tuple(sample_id, markdup.bam)
 * Output: tuple(sample_id, flagstat.txt)
 */
process SAMTOOLS_FLAGSTAT {
    tag "${sample_id}"

    publishDir "${params.outdir}/markdup_bams", mode: 'copy', pattern: "*.flagstat.txt"

    input:
    tuple val(sample_id), path(markdup_bam)

    output:
    tuple val(sample_id), path("${sample_id}.flagstat.txt")

    script:
    """
    samtools flagstat "${markdup_bam}" > ${sample_id}.flagstat.txt
    """
}

// process FLAGSTAT_SUMMARY_MATRIX {
//     tag "flagstat_summary"

//     publishDir "${params.outdir}/markdup_bams", mode: 'copy', pattern: "samtools_flagstat_matrix.tsv"

//     input:
//     val flagstat_lines // a single value: list or string, depending on how you use it

//     output:
//     path "samtools_flagstat_matrix.tsv"

//     script:
//     """
//     python /Users/massal01/0_base_scr/cytoDRIP_align/cytoDRIP_HM/report_code/agg_flagstat_matrix.py ${flagstat_lines} samtools_flagstat_matrix.tsv
//     """
// }