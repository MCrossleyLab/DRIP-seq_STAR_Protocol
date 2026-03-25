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

        def fastq1 = file(m.fastq1.toString())
        def fastq2 = m.fastq2 ? file(m.fastq2.toString()) : null

        // species parsed into a list
        def raw_species = (m.species ?: '').toString().trim()
        def species_list = raw_species ? raw_species.split(/[;,]/)*.trim().findAll { it } : []

        // emit tuple: mandatory fields + full metadata map
        tuple( sample_id, fastq1, fastq2, species_list, m )
    }

    // aligned_input_ch = samples_ch.map { sid, fastq1, fastq2, species, treatment, replicate, tissue, condition -> tuple(sid, fastq1, read2, species) }
    aligned_input_ch = samples_ch.map { sid, fastq1, fastq2, species, m -> tuple(sid, fastq1, fastq2, species) }
    
    // 2) Bowtie reads
    aligned_ch = bowtie_reads(aligned_input_ch) // (sid, fastq1, fastq2, species) >> (sid, species, aligned_bam)
    
    // 3) Mark duplicates
    markdup_out = markdup_bams(aligned_ch) // (sid, species, aligned_bam) >> (sid, species, markdup_bam)

    // 4) Get metadata
    meta_ch = samples_ch.map{ sid, fastq1, fastq2, species, m -> tuple(sid, m.treatment, m.replicate) }
    
    // 4a) combine metadata
    bam_with_meta_ch = markdup_out.markdup_ch
        .combine(meta_ch, by: 0) // tuple(sid, species, markdup_bam, treatment, replicate)
    
    // 4b) make bigWigs
    bam_norm_out = bam_to_bigWig(bam_with_meta_ch)

}

