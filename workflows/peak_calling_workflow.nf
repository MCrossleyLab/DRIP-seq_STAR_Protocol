#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { MACS2_CALLPEAK } from '../modules/peak_calling_modules.nf'


workflow peak_calling {
    take:
    bam_with_meta_ch

    main:
    bam_records = bam_with_meta_ch.toList()
    explicit_requests_ch = bam_records
        .map { recs ->
            def bam_by_sid = recs.collectEntries { rec ->
                [(rec[0]): rec[2]]
            }

            params.contrasts.collect { contrast_name, sides ->
                def ip_bams = sides[0].collect { sid -> bam_by_sid[sid] }
                def in_bams = sides[1].collect { sid -> bam_by_sid[sid] }
                
                [contrast_name, ip_bams, in_bams]
            }
        }
        .flatMap { it }
    
    callpeaks_ch = MACS2_CALLPEAK(explicit_requests_ch)

    emit:
    callpeaks_ch.broadpeaks
}