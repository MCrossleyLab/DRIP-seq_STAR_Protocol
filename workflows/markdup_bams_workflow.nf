#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { SAMTOOLS_COLLATE;
          SAMTOOLS_FIXMATE;
          SAMTOOLS_SORT;
          SAMTOOLS_FLAGSTAT;
          SAMTOOLS_MARKDUP_INDEX } from '../modules/markdup_bams_modules.nf'

workflow markdup_bams {
    take:
    aligned_ch

    main:
    aligned_input_ch    = aligned_ch.map { sid, species_list, bam -> tuple(sid, bam) }
    collate_ch          = SAMTOOLS_COLLATE(aligned_input_ch)
    fixmate_ch          = SAMTOOLS_FIXMATE(collate_ch)      
    sorted_ch           = SAMTOOLS_SORT(fixmate_ch)         
    markdup_no_sp_ch    = SAMTOOLS_MARKDUP_INDEX(sorted_ch) 
    
    species_meta_ch = aligned_ch.map { sid, species_list, bam -> tuple(sid, species_list) }
    markdup_ch = species_meta_ch.join(markdup_no_sp_ch)
    
    flagstat_ch = SAMTOOLS_FLAGSTAT(markdup_no_sp_ch)  
    flagstat_list_ch = flagstat_ch
        .map { sid, f -> "${sid}\t${f}" }
        .toSortedList { a, b ->
            def runA = a.split('\t')[0]
            def runB = b.split('\t')[0]
            runA <=> runB
        }
        .collect()
        .map { lines -> lines.join(' ') }

    emit:
    markdup_ch
}
