
#DRIP-seq paired end bioinformatics analysis

#Basic alignment to genome and making genome browser tracks

#Quality control check with fastqc
find *.fastq.gz | grep -v 'Und' | xargs fastqc -o fastqc/ --threads 10

#Trim Tru-seq illumina adapters
cutadapt -a AGATCGGAAGA -A AGATCGGAAGA -m 40 -o out.1.trimmed.fastq -p out.2.trimmed.fastq reads.1.fastq reads.2.fastq

# Building bowtie index
bowtie2 index mm10.fa mm10

#Alignment to mm10 genome using bowtie2 end to end settings
bowtie2 -p 10 --no-unal --no-mixed --no-discordant -x mm_10 -1 reads1.fastq -2 reads2.fastq | samtools view -bSh > alignment.bam

#Sort and index bam files
ml load samtools
samtools sort -@ 8 -o alignment.sorted.bam alignment.bam
samtools index alignment.sorted.bam

# Remove duplicates, optional as duplicates are removed with gsort -u below, use these bam files for peak calling
java -jar picard.jar MarkDuplicates REMOVE_DUPLICATES=true I=alignment.sorted.bam O=alignment.sorted.rmdup.bam M=metrics.txt

#Since it's paired end data, convert from BAM to BEDPE
cat files_names.txt | sed 's/\rmdup.bam//' | xargs -I {} -P10 sh -c 'bedtools bamtobed -bedpe -i {}.bam > {}.bedpe'

#Convert fragments to bed files, take columns 1 (chr R1), 2 (start R1), and 6 (end R2)
find *.bedpe | sed 's/\.bedpe//' | xargs -I {} -c 'cut -f 1,2,6 {}.bedpe > {}.bed'

#Sort and remove any duplicates (note, picard step above not strictly required)
find *.bed | sed 's/\.bed//' | xargs -I {} -P10 sh -c 'gsort --parallel=8 -u -k1,1 -k2,2n -k3,3n {}.bed > {}.sorted.bed'

#Convert to bedGraph
find *.sorted.bed | sed 's/\.sorted.bed//' | xargs -I {} -P10 sh -c 'bedtools genomecov -bg -g mm10.chrom.sizes -i {}.sorted.bed > {}.bedGraph'
#Count the reads for normalization, convert to RPM
wc -l *.bed > read_count.txt

9069956	V1_Input-1_Hepa.sorted.bed	9.069956
8979683	V1_Input-1_L.sorted.bed	8.979683
7755421	V1_Input-2_Hepa.sorted.bed	7.755421
11282307	V1_Input-2_L.sorted.bed	11.282307
14300068	V1_RNH--1_Hepa.sorted.bed	14.300068
14098233	V1_RNH--1_L.sorted.bed	14.098233
26153489	V1_RNH--2_Hepa.sorted.bed	26.153489
26129303	V1_RNH--2_L.sorted.bed	26.129303
12217100	V1_RNHplus-1_Hepa.sorted.bed	12.2171
8979683	V1_RNHplus-1_L.sorted.bed	8.979683
20554130	V1_RNHplus-2_Hepa.sorted.bed	20.55413
12257706	V1_RNHplus-2_L.sorted.bed	12.257706

#Normalize bedGraphs to RPM for each sample file
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4/9.069956)}' V1_Input-1_Hepa.bedGraph > V1_Input-1_Hepa.bedGraph.norm
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4/8.979683)}' V1_Input-1_L.bedGraph > V1_Input-1_L.bedGraph.norm
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4/7.755421)}' V1_Input-2_Hepa.bedGraph > V1_Input-2_Hepa.bedGraph.norm
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4/11.282307)}' V1_Input-2_L.bedGraph > V1_Input-2_L.bedGraph.norm

awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4/14.300068)}' V1_RNH--1_Hepa.bedGraph > V1_RNH--1_Hepa.bedGraph.norm
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4/14.098233)}' V1_RNH--1_L.bedGraph > V1_RNH--1_L.bedGraph.norm
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4/26.153489)}' V1_RNH--2_Hepa.bedGraph > V1_RNH--2_Hepa.bedGraph.norm
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4/26.129303)}' V1_RNH--2_L.bedGraph > V1_RNH--2_L.bedGraph.norm

awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4/12.2171)}' V1_RNHplus-1_Hepa.bedGraph > V1_RNHplus-1_Hepa.bedGraph.norm
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4/8.979683)}' V1_RNHplus-1_L.bedGraph > V1_RNHplus-1_L.bedGraph.norm
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4/20.55413)}' V1_RNHplus-2_Hepa.bedGraph > V1_RNHplus-2_Hepa.bedGraph.norm
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4/12.257706)}' V1_RNHplus-2_L.bedGraph > V1_RNHplus-2_L.bedGraph.norm

#Convert to bigwig files
find *.norm | sed 's/\.bedGraph.norm//' | xargs -I {} -P 10 /util/bedGraphToBigWig {}.bedGraph.norm ../mm10.chrom.sizes {}.bw


#Alternative normalization straight from the bam files
bamCoverage --bam V1_Input-2.alignment_mark_duplicates.bam -o V1_Input-2.SeqDepthNorm.bedgraph --outFileFormat bedgraph --binSize 1 --numberOfProcessors 40 --normalizeUsing CPM
sort -k1,1 -k2,2n -k3,3n V1_Input-2.SeqDepthNorm.bedgraph > V1_Input-2.SeqDepthNorm.sorted.bedgraph
ml load ucsc-tools
bedGraphToBigWig V1_Input-2.SeqDepthNorm.sorted.bedgraph mm10.chrom.sizes V1_Input-2.bigwig


#Peak calling and peak-centered analysis
#Merge biological replicates for input and IP samples before peak calling

samtools merge merged/liver_IN_merged.bam V1_Input-1_L.rmdup.bam V1_Input-2_L.alignment.rmdup.bam
samtools merge merged/liver_IP_merged.bam V1_RNH--1_L_alignment.rmdup.bam V1_RNH--2_L.alignment.sorted.rmdup.bam
samtools merge merged/hepa_IN_merged.bam V1_Input-2_Hepa_alignment.sorted.rmdup.bam V1_Input-1_Hepa.alignment.sorted.rmdup.bam
samtools merge merged/hepa_IP_merged.bam V1_RNH--1_Hepa_alignment.sorted.rmdup.bam V1_RNH--2_Hepa.alignment.rmdup.bam

#Call peaks against input using MACS2
macs2 callpeak --broad -g mm -f BAMPE -c merged/liver_IN_merged.bam -t merged/liver_IP_merged.bam --outdir ./liver_peak_calls -n liver_merged_IP_vs_IN 2> ./liver_merged_IP_vs_IN_broadpeaks.log
macs2 callpeak --broad -g mm -f BAMPE -c merged/hepa_IN_merged.bam -t merged/hepa_IP_merged.bam --outdir ./hepa_peak_calls -n hepa_merged_IP_vs_IN 2> ./hepa_merged_IP_vs_IN_broadpeaks.log

#convert to bed files
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4,$5,$6)}' hepa_merged_IP_vs_IN_peaks.broadPeak > hepa_peaks.bed
awk 'BEGIN{OFS="\t"}{print($1,$2,$3,$4,$5,$6)}' liver_merged_IP_vs_IN_peaks.broadPeak > liver_peaks.bed

#About 28,000 peaks in hepa
#About 14,000 in liver
#Load peak bed files on genome browser like IGV with bigwig files to check by insepction that they look reasonable 

#DRIP peak annotation
#Now overlap peaks with mm10 genes to get those that are + or - stranded and also genic vs intergenic

#First do a stranded slop on full mm10 genes file - to include 3kb upstream of TSS and 3Kb downstream of TES. And define these as genic regions
#Download from UCSC genome browser
cut -f 1,2,3,4,5,6 mm10_gencodev25.basic.bed > mm10_v25_basic.6col.bed
bedtools slop -i mm10_v25_basic.6col.bed -g mm10.chrom.sizes -b 3000 > mm10_v25_3kb_slop_genes.bed

#Intersect mm10 genes with DRIP peaks
bedtools intersect -wo -b hepa_annotated_peaks.bed  -a mm10_v25_3kb_slop_genes.bed > hepa_intersect_v25_genes.tsv
bedtools intersect -wo -b liver_annotated_peaks.bed  -a mm10_v25_3kb_slop_genes.bed > liver_intersect_v25_genes.tsv

#assign gene strand to peak strand

#Use deeptools to generate plots over TSS and over peak centers, binsize 75
#DRIP peak centers
cat sample_names.txt | xargs -I % computeMatrix reference-point -R hepa_annotated_peaks.bed -S bigwigs_re/%.bw --outFileName deeptools_hepa_bin75/%.matrix -bs 75 -a 7500 -b 7500 -p 8 --referencePoint center
cat sample_names.txt | xargs -I % computeMatrix reference-point -R liver_annotated_peaks.bed -S bigwigs_re/%.bw --outFileName deeptools_liver_bin75/%.matrix -bs 75 -a 7500 -b 7500 -p 8 --referencePoint center

#Over mm10 genes TSS and TES
cat sample_names.txt | xargs -I % computeMatrix reference-point -R mm10_v25_basic.6col.bed -S bigwigs_re/%.bw --outFileName deeptools_tss_bin75/%.matrix -bs 75 -a 7500 -b 7500 -p 8 --referencePoint TSS
cat sample_names.txt | xargs -I % computeMatrix reference-point -R mm10_v25_basic.6col.bed -S bigwigs_re/%.bw --outFileName deeptools_tes_bin75/%.matrix -bs 75 -a 7500 -b 7500 -p 8 --referencePoint TES


# Nucleotide content over peaks
bedtools nuc -fi mm10.fa -bed hepa_annotated_peaks.bed   > hepa_annotated_peaks.bed.nuc_content.txt
bedtools nuc -fi mm10.fa -bed liver_annotated_peaks.bed  > liver_annotated_peaks.bed.nuc_content.txt


#Non-B DNA structures at DRIP peaks

#https://ncifrederick.cancer.gov/bids/ftp/?nonb

#Have a database of predicted non-B-form DNA structures
#download gz file, double click in finder to extract and also extract tarball
#Get folder with type of secondary structure separated for each chromosomes
#-All are gff or tsv files (12 col) with chromosome for each structure type
#Use TSV files with following header
#Sequence_name	Source	Type	Start	Stop	Length	Score	Strand	Repeat	Spacer	Permutations	Subset	Composition	Sequence
#chr10	ABCC	Mirror_Repeat	3001161	3001183	23	NA	+	11	1	1	1	9A/1C/0G/1T	taaaaaacaaacaaacaaaaaat
#Structure classes are:
#APR - a phased repeat
#DR - direct repeat
#GQ - G-quad
#IR - inverted repeat
#MR - mirror repeat
#STR - short tandem repeat
#Z - Z DNA
#Cat files from each chromosome for each structure type:
#remove chrMT for each type (only a couple, but may skew things)
#e.g.
rm chrMT_Z.tsv 
#Then combine

cat *Z* | grep -v 'Sequence_name' > Z.tsv
cat *APR* | grep -v 'Sequence_name' > APR.tsv
cat *GQ* | grep -v 'Sequence_name' > GQ.tsv
cat *IR* | grep -v 'Sequence_name' > IR.tsv
cat *MR* | grep -v 'Sequence_name' > MR.tsv
cat *STR* | grep -v 'Sequence_name' > STR.tsv
cat *DR* | grep -v 'Sequence_name' > DR.tsv


#1)Get bed files of each structure:
#2)Then add column with 1, this marks the presence of the nonB DNA structure
#3)Convert to bigwig
#4)Deeptools gene metaplot and around peak centers (try peak centers 1st)

#1) 
cat nonbdna_types.txt | xargs -I {} -P6 sh -c 'cut -f 1,4,5 {}.tsv | sort -k1,1 -k2,2n -k3,3n > {}.bed'

#2)
cat APR.bed | awk '{ print $1 "\t" $2 "\t" $3 "\t1"}' > APR.4col.bed
cat DR.bed | awk '{ print $1 "\t" $2 "\t" $3 "\t1"}' > DR.4col.bed
cat GQ.bed | awk '{ print $1 "\t" $2 "\t" $3 "\t1"}' > GQ.4col.bed
cat IR.bed | awk '{ print $1 "\t" $2 "\t" $3 "\t1"}' > IR.4col.bed
cat MR.bed | awk '{ print $1 "\t" $2 "\t" $3 "\t1"}' > MR.4col.bed
cat STR.bed | awk '{ print $1 "\t" $2 "\t" $3 "\t1"}' > STR.4col.bed
cat Z.bed | awk '{ print $1 "\t" $2 "\t" $3 "\t1"}' > Z.4col.bed

#Download mm9 chrom.sizes

wget --timestamping \
        'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes' \
        -O mm9.chrom.sizes

#3) 

ls -1 *4col.bed | sed 's/.4col.bed//' | xargs -I {} -P6 sh -c 'bedtools genomecov -bg -g mm9.chrom.sizes -i {}.4col.bed > {}.bedGraph'
ls -1 *4col.bed | sed 's/.4col.bed//' | xargs -I {} -P6 sh -c  'util/bedGraphToBigWig {}.bedGraph mm9.chrom.sizes {}.bw'

#4) 
#Get control DRIP regions - in same genes as DRIP peaks (i.e. expression matched)

#-Get DRIP peaks and slop 3000 to make sure any other peaks will be at least 3kb away
#-Then intersect these with genes that contain a DRIP peak
#-Then do random shuffle from within these genes excluding the slopped peaks
#-Use LiftOver to get these to mm9 (repeat bigwigs are in mm9)

bedtools slop -i hepa_peaks.mm9.bed -g mm9.chrom.sizes -b 3000 > hepa_peaks_slop_3kb.mm9.bed
bedtools slop -i liver_peaks.mm9.bed -g mm9.chrom.sizes -b 3000 > liver_peaks_slop_3kb.mm9.bed

bedtools slop -i 3T3_DRIP_peaks.mm9.bed -g mm9.chrom.sizes -b 3000 > t3_peaks_slop_3kb.mm9.bed
bedtools slop -i E14_DRIP_peaks.mm9.bed -g mm9.chrom.sizes -b 3000 > e14_peaks_slop_3kb.mm9.bed

#where -g is mm9 chrom sizes
#Need to download:

wget 'ftp://hgdownload.cse.ucsc.edu/goldenPath/mm9/bigZips/mm9.chrom.sizes' -O mm9.chrom.sizes

bedtools shuffle -incl mm9_genes_6col.bed -excl hepa_peaks_slop_3kb.mm9.bed -i hepa_peaks.mm9.bed -g mm9.chrom.sizes > hepa_peaks_shuffle.1.bed
bedtools shuffle -incl mm9_genes_6col.bed -excl liver_peaks_slop_3kb.mm9.bed -i liver_peaks.mm9.bed -g mm9.chrom.sizes > liver_peaks_shuffle.1.bed

#5) Now run in deeptools to get repeat density around DRIP peaks and shuffled controls

cat sample_names2.txt | xargs -I % computeMatrix reference-point -R liver_annotated_peaks.bed -S bigwigs_re/%.bw --outFileName deeptools/liver_center_re/%.matrix -bs 10 -a 5000 -b 5000 -p 8 --referencePoint center


