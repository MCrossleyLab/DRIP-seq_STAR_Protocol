#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { COUNT_READS_RPM;
          BAMCOVERAGE } from '../modules/bam_to_bigWig_modules.nf'


workflow bam_to_bigWig {
    take:
    bam_ch

    main:
    counts_ch = bam_ch | COUNT_READS_RPM

    per_sample_scale_ch = counts_ch
        .groupTuple(by: [0,1,2,3])
        .map { row ->
            def (sid, species, treat, rep, non_dup_list) = row
            def non_dup_mean = non_dup_list.collect { it as long }.sum() / non_dup_list.size()

            def rpm = non_dup_mean ? (non_dup_mean / 1_000_000.0) : 0.0
            def rpm_scale = rpm ? (1.0 / rpm) : 0.0

            tuple(sid, species, treat, rep, non_dup_mean, rpm_scale)
        }

    scale_key_ch = bam_ch
        .map { sid, species, bam_file, treat, rep ->
                tuple(sid, species, treat, rep, bam_file)
            }
        .combine(per_sample_scale_ch, by: [0,1,2,3])
        .flatMap { sid, species_list, treat, rep, bam_file, non_dup_mean, rpm_scale ->
            if (!species_list || species_list.isEmpty()) {
                [ tuple(sid, 'nospecies', bam_file, treat, rep, non_dup_mean, rpm_scale) ]
            } else {
                species_list.collect { sp -> tuple(sid, sp, bam_file, treat, rep, non_dup_mean, rpm_scale) }
            }
        }
    
    bamCov_bw_ch = scale_key_ch | BAMCOVERAGE

    emit:
    bamCov_bw_ch
}
