#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { COUNT_READS_RPM;
          WRITE_BAM_RPM_TABLE; 
          BAMCOVERAGE } from '../modules/bam_to_bigWig_modules.nf'


workflow bam_to_bigWig {
    take:
    bam_ch // tuple(sid, species_list, bam_file, treatment, replicate)

    main:
    //------------------------------------------------------------------
    // 1) non-duplicate counts
    //------------------------------------------------------------------
    counts_ch = bam_ch | COUNT_READS_RPM
    // (sid, species, treat, rep, non_dup)
    
    // counts_ch | view
    // [SRR28361292, [mm10], HEPA1-6, rep3, 56685144]
    //------------------------------------------------------------------
    // 2) mean non-duplicate reads per sample/species/treatment/replicate
    //------------------------------------------------------------------
    per_sample_scale_ch = counts_ch
        // .map { sid, species, treat, rep, non_dup ->
        //     tuple(sid, species, treat, rep, non_dup)
        // }
        .groupTuple(by: [0,1,2,3])   // group by (sid,species,treat,rep)
        .map { row ->
            def (sid, species, treat, rep, non_dup_list) = row
            def non_dup_mean = non_dup_list.collect { it as long }.sum() / non_dup_list.size()

            def rpm = non_dup_mean ? (non_dup_mean / 1_000_000.0) : 0.0
            def rpm_scale = rpm ? (1.0 / rpm) : 0.0

            tuple(sid, species, treat, rep, non_dup_mean, rpm_scale)
        }
    // per_sample_scale_ch | view 
    // [SRR28361297, [mm10], Liver, rep6, 26322062, 0.0379909446]
    //------------------------------------------------------------------
    // 3) Summary table (only non_dup + factor)
    //------------------------------------------------------------------
    bam_rpm_table_ch = per_sample_scale_ch
        .map { sid, species, treat, rep, non_dup_mean, rpm_scale ->
            "${sid}\t${species}\t${treat}\t${rep}\t${non_dup_mean}\t${rpm_scale}"
        }
        .toSortedList { a, b ->
            def aID = a.split('\\t')[0]
            def bID = b.split('\\t')[0]
            aID <=> bID
        }
        .collect()

    bam_rpm_table_ch | WRITE_BAM_RPM_TABLE

    // //------------------------------------------------------------------
    // // 4) Join per-sample non_dup and RPM factor back to each strand bedGraph
    // //------------------------------------------------------------------
    // rpm_key_ch = per_sample_counts_ch
    //     .map { sid, species, treat, rep, non_dup_mean, rpm_scale ->
    //         tuple(sid, species, treat, rep, non_dup_mean, rpm_scale)
    //     }

    // bedgraph_key_ch = bedgraph_with_meta_ch
    //     .map { sid, species, strand, bedgraph, treat, rep ->
    //         tuple(sid, species, treat, rep, strand, bedgraph)
    //     }

    // bedgraph_with_rpm_ch = bedgraph_key_ch
    //     .combine(rpm_key_ch, by: [0,1,2,3])
    //     .map { sid, species, treat, rep, strand, bedgraph, non_dup_mean, rpm_scale ->
    //         tuple(sid, species, strand, bedgraph, treat, rep, non_dup_mean, rpm_scale)
    //     }
    // // bedgraph_with_rpm_ch | view 
    
    // //------------------------------------------------------------------
    // // 5) Apply RPM normalization (and keep non_dup in env)
    // //------------------------------------------------------------------
    // norm_bam_rpm_ch = bedgraph_with_rpm_ch | NORM_BEDGRAPH_RPM_FROM_BAM
    
    //------------------------------------------------------------------
    // 6) Direct bigWig from BAM via deepTools (per strand)
    //------------------------------------------------------------------
    scale_key_ch = bam_ch //tuple(sid, species, bam_file, treatment, replicate)
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

    // scale_key_ch | view 
    
    bamCov_bw_ch = scale_key_ch | BAMCOVERAGE

    emit:
    // strand_counts_ch
    // per_sample_counts_ch
    // bam_rpm_table_ch
    // norm_bam_rpm_ch
    bamCov_bw_ch
}
