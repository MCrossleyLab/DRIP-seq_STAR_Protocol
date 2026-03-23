#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { bowtie_reads        } from './workflows/bowtie_reads_workflow.nf'
include { markdup_bams        } from './workflows/markdup_bams_workflow.nf'
include { bam_to_bigWig       } from './workflows/bam_to_bigWig_workflow.nf'


workflow {

    // 1) Read samplesheet -> channel of tuples
    ch_meta = Channel
    .fromPath(params.samplesheet)
    .splitCsv(header: true)   // row is a Map: [col1: val1, col2: val2, ...]
    .map { row ->
        // normalise keys if you like
        def m = [:]
        row.each { k, v -> m[k.toString().toLowerCase()] = v }
        m
    }
    
    samples_ch = ch_meta.map { m ->

        // required fields
        def sample_id = m.sample_id.toString()

        def read1 = file(m.read1.toString())
        def read2 = m.read2 ? file(m.read2.toString()) : null

        // species parsed into a list
        def raw_species = (m.species ?: '').toString().trim()
        def species_list = raw_species ? raw_species.split(/[;,]/)*.trim().findAll { it } : []

        // emit tuple: mandatory fields + full metadata map
        tuple( sample_id, read1, read2, species_list, m )
    }

    // samples_ch | view

    // aligned_input_ch = samples_ch.map { sid, read1, read2, species, treatment, replicate, tissue, condition -> tuple(sid, read1, read2, species) }
    aligned_input_ch = samples_ch.map { sid, read1, read2, species, m -> tuple(sid, read1, read2, species) }
    
    // 2) Bowtie reads
    aligned_ch = bowtie_reads(aligned_input_ch) // (sid, read1, read2, species) >> (sid, species, aligned_bam)
    
    // 3) Mark duplicates
    markdup_out = markdup_bams(aligned_ch) // (sid, species, aligned_bam) >> (sid, species, markdup_bam)

    // 6a) norm from bedGraph
    // meta_ch = samples_ch.map{ sid, read1, read2, species, treatment, replicate, tissue, condition -> tuple(sid, treatment, replicate) }
    meta_ch = samples_ch.map{ sid, read1, read2, species, m -> tuple(sid, m.treatment, m.replicate) }
    
    // 6b) norm from bams
    bam_with_meta_ch = markdup_out.markdup_ch
        .combine(meta_ch, by: 0) // tuple(sid, species, markdup_bam, treatment, replicate)
    
    // Run RPM-from-BAM workflow
    bam_norm_out = bam_to_bigWig(bam_with_meta_ch)
    

    // 3) Join BAMs with metadata by sample_id
    design_mat_ch = samples_ch.map { sample_id, read1, read2, spec, m ->
        tuple(
            id:   sample_id,
            [
                species:   spec,    // mm10
                treatment: m.treatment.toString(),  // input/ip/ip.RNAseH1
                rep:       m.replicate.toString(),
                tissue:    m.tissue.toString(),     // Liver / HEPA1-6
                cond:      m.condition.toString()
            ]
        )
    }
    
    bam_ch = markdup_out.markdup_ch.map { sid, spec, bam ->
        def species = spec[0]    // mm10
        tuple(
            id: sid,
            [
                species: spec,
                bam:     bam
            ]
        )
    }

    metadata_ch = bam_ch
        .combine(design_mat_ch, by: 0)
        .map { lst -> lst.inject([:]) { acc, m -> acc + m } }  // flatten dicts into one map

    //  (meta, [bam1, bam2, ...]) grouped per (species, tissue, treatment)
    grouped_ch = metadata_ch
    .map { m ->
        // group key as small map or tuple
        def meta = [
            species: m.species,
            tissue:  m.tissue,
            tret:    m.treatment
        ]
        tuple(meta, m.bam)
    }
    .groupTuple()

    // grouped_ch: (meta, bams)

    liver_input_ch = grouped_ch.filter { meta, bams ->
        meta.tissue == 'Liver' && meta.tret == 'input'
    }

    liver_ip_ch = grouped_ch.filter { meta, bams ->
        meta.tissue == 'Liver' && meta.tret == 'ip'
    }

    liver_rnaseh1_ch = grouped_ch.filter { meta, bams ->
        meta.tissue == 'Liver' && meta.tret == 'ip.RNAseH1'
    }

    hepa_input_ch = grouped_ch.filter { meta, bams ->
        meta.tissue == 'HEPA1-6' && meta.tret == 'input'
    }

    hepa_ip_ch = grouped_ch.filter { meta, bams ->
        meta.tissue == 'HEPA1-6' && meta.tret == 'ip'
    }

    hepa_rnaseh1_ch = grouped_ch.filter { meta, bams ->
        meta.tissue == 'HEPA1-6' && meta.tret == 'ip.RNAseH1'
    }


    // branches.liver_ip
    // .map { meta, bams -> tuple(meta.species, meta.tissue, meta.tret, bams.flatten()) }
    // | view

    
    // // TODO: TEST THIS
    // // #Alternative normalization straight from the bam files
    // // bamCoverage --bam V1_Input-2.alignment_mark_duplicates.bam -o V1_Input-2.SeqDepthNorm.bedgraph --outFileFormat bedgraph --binSize 1 --numberOfProcessors 40 --normalizeUsing CPM
    // // sort -k1,1 -k2,2n -k3,3n V1_Input-2.SeqDepthNorm.bedgraph > V1_Input-2.SeqDepthNorm.sorted.bedgraph
    // // ml load ucsc-tools
    // // bedGraphToBigWig V1_Input-2.SeqDepthNorm.sorted.bedgraph mm10.chrom.sizes V1_Input-2.bigwig
    
    // // Peak calling
    // #Peak calling and peak-centered analysis
    // #Merge biological replicates for input and IP samples before peak calling

    // samtools merge merged/liver_IN_merged.bam V1_Input-1_L.rmdup.bam V1_Input-2_L.alignment.rmdup.bam
    // samtools merge merged/liver_IP_merged.bam V1_RNH--1_L_alignment.rmdup.bam V1_RNH--2_L.alignment.sorted.rmdup.bam
    // samtools merge merged/hepa_IN_merged.bam V1_Input-2_Hepa_alignment.sorted.rmdup.bam V1_Input-1_Hepa.alignment.sorted.rmdup.bam
    // samtools merge merged/hepa_IP_merged.bam V1_RNH--1_Hepa_alignment.sorted.rmdup.bam V1_RNH--2_Hepa.alignment.rmdup.bam

    // #Call peaks against input using MACS2
    // macs2 callpeak --broad -g mm -f BAMPE -c merged/liver_IN_merged.bam -t merged/liver_IP_merged.bam --outdir ./liver_peak_calls -n liver_merged_IP_vs_IN 2> ./liver_merged_IP_vs_IN_broadpeaks.log
    // macs2 callpeak --broad -g mm -f BAMPE -c merged/hepa_IN_merged.bam -t merged/hepa_IP_merged.bam --outdir ./hepa_peak_calls -n hepa_merged_IP_vs_IN 2> ./hepa_merged_IP_vs_IN_broadpeaks.log

    // #convert to bed files
    // awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4,$5,$6)}' hepa_merged_IP_vs_IN_peaks.broadPeak > hepa_peaks.bed
    // awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4,$5,$6)}' liver_merged_IP_vs_IN_peaks.broadPeak > liver_peaks.bed
}

