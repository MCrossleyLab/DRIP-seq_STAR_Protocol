# DRIP‑seq pipeline

This repository contains the analysis code used for DRIP‑seq experiments genereated by Puzzo et al. [https://www.sciencedirect.com/science/article/pii/S1525001624006609?via%3Dihub] and described in the STAR Protocol manuscript.

## inputs

* samplesheet.csv   : headings (sample_id,read1,read2,species,treatment,replicate)
* parameters.yml    : pipeline parameters, such as: output_dir, adapter_fwd, adapter_fwd, ...

## worflow

1. Generate aligned BAMs from fastq
2. Calculate normalisation factors
3. generate BigWig track files

## pipeline requirments

```sh
python=3.11.11
fastqc=0.12.1
cutadapt==5.2
bowtie2=2.5.4
samtools=1.6
pybedtools=0.12.0
bedtools=2.31.1
deeptools=3.5.6
```
# DRIP-seq_STAR_Protocol
