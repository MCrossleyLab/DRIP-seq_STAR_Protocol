#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include { CUTADAPT_ADAPTER_TRIM; 
          BOWTIE2_ALIGN } from '../modules/bowtie_reads_modules.nf'

def findIndexFromSpecies(List species_list, String idx_dir, String suffix = "") {
    def uniq = species_list as List
    uniq = uniq.collect { it?.toString()?.trim() }.findAll { it }      // non-empty strings
    uniq = uniq.unique()

    if (!uniq) {
        throw new IllegalArgumentException("Empty species_list for sample")
    }

    def keys = []

    if (uniq.size() == 1) {
        keys << uniq[0]
    } else {
        for (int i = 0; i < uniq.size(); i++) {
            for (int j = i + 1; j < uniq.size(); j++) {
                def a = uniq[i]
                def b = uniq[j]
                keys << "${a}_${b}"
                keys << "${b}_${a}"
            }
        }
        keys << uniq.join('_')
        keys << uniq.sort().join('_')
    }

    keys = keys.unique()

    for (key in keys) {
        def basename = suffix ? "${key}_${suffix}" : key
        def idx_base = "${idx_dir}/${basename}"

        def dir = new File(idx_dir)
        def bt2_files = dir.listFiles()?.findAll { f ->
            f.name.startsWith("${basename}.") && f.name.endsWith(".bt2")
        }

        if (bt2_files && !bt2_files.isEmpty()) {
            return idx_base
        }
    }

    throw new IllegalArgumentException(
        "No Bowtie2 index found in ${idx_dir} for species_list=${uniq} with keys=${keys}"
    )
}


workflow bowtie_reads {

    take:
    samples_ch

    main:
    raw_pairs_ch     = samples_ch.map { sid, r1, r2, sp -> tuple(sid, file(r1), file(r2)) }
    raw_pairs_ch | CUTADAPT_ADAPTER_TRIM
    
    species_ch = samples_ch.map { sid, r1, r2, sp -> tuple(sid, sp) }
    joined_ch = CUTADAPT_ADAPTER_TRIM.out.sid_fastq_files.join(species_ch)
    
    bowtie_input_ch = joined_ch.map { sid, r1, r2, species_list ->
        def idx_base = findIndexFromSpecies(
            species_list as List,
            params.bowtie2_index_dir as String,
            params.bowtie2_index_suffix as String
        )
        tuple(sid, r1, r2, params.bowtie2_extra_params, idx_base, species_list)
    }
    
    BOWTIE2_ALIGN(bowtie_input_ch)
    
    emit:
    BOWTIE2_ALIGN.out.sid_bam_file
}
