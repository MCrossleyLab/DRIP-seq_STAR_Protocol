#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { SAMTOOLS_COLLATE;
          SAMTOOLS_FIXMATE;
          SAMTOOLS_SORT;
          SAMTOOLS_FLAGSTAT;
          SAMTOOLS_MARKDUP_INDEX } from '../modules/markdup_bams_modules.nf'

workflow markdup_bams {
    take:
    aligned_ch    // (sample_id, species_list, bam)

    main:
    aligned_input_ch    = aligned_ch.map { sid, species_list, bam -> tuple(sid, bam) }
    collate_ch          = SAMTOOLS_COLLATE(aligned_input_ch)    // (sample_id, bam)
    fixmate_ch          = SAMTOOLS_FIXMATE(collate_ch)          // (sample_id, bam)
    sorted_ch           = SAMTOOLS_SORT(fixmate_ch)             // (sample_id, bam)
    markdup_no_sp_ch    = SAMTOOLS_MARKDUP_INDEX(sorted_ch)     // (sample_id, bam)
    
    // restore species info
    species_meta_ch = aligned_ch.map { sid, species_list, bam -> tuple(sid, species_list) }
    markdup_ch = species_meta_ch.join(markdup_no_sp_ch) // (sample_id, species_list, markup_bam)
    
    flagstat_ch = SAMTOOLS_FLAGSTAT(markdup_no_sp_ch)         // input (sample_id, bam) => emits (sample_id, flagstat_file)
    flagstat_list_ch = flagstat_ch
        .map { sid, f -> "${sid}\t${f}" }
        .toSortedList { a, b ->                      // sort the strings
            def runA = a.split('\t')[0]
            def runB = b.split('\t')[0]
            runA <=> runB                             // lexicographic by runID
        }
        .collect()
        .map { lines -> lines.join(' ') }

    // flagstat_matrix_ch = FLAGSTAT_SUMMARY_MATRIX(flagstat_list_ch)

    emit:
    markdup_ch          // (sample_id, species_list, markup_bam)
    // flagstat_matrix_ch
}
