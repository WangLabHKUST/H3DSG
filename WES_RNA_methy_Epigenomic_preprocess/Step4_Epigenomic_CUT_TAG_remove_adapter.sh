#!/bin/bash
#SBATCH --mail-user=shanghaixia@ust.hk #Update your email address
#SBATCH --mail-type=end
#SBATCH -p cpu-share
#SBATCH -N 1 -n 20
dir=/scratch/PI/jgwang/haixia/projects/SP_CUT_Tag
ID=${1}
cd /scratch/PI/jgwang/haixia/projects/SP_CUT_TAG_new_mapping
#mkdir ${ID}_cut_tag_remove_adapter
cd ${ID}_cut_tag_remove_adapter
tar=/scratch/PI/jgwang/haixia/projects/SP_CUT_TAG_new_mapping/${ID}_cut_tag_remove_adapter
Ref1=/scratch/PI/jgwang/haixia/projects/SP_Chipseq/bowtie2_indexes/hg19

java=/home/shanghaixia/miniconda3/bin/java
samtools=/home/shanghaixia/miniconda3/bin/samtools
bowtie2=/home/shanghaixia/miniconda3/bin/bowtie2
macs2=/home/shanghaixia/miniconda3/bin/macs2
bioawk=/home/shanghaixia/miniconda3/bin/bioawk
bedtools=/home/shanghaixia/miniconda3/bin/bedtools
sort=/usr/bin/sort
bedClip=/home/shanghaixia/miniconda3/bin/bedClip
bedGraphToBigWig=/home/shanghaixia/miniconda3/bin/bedGraphToBigWig
fastp=/home/shanghaixia/miniconda3/bin/fastp
fastqc=/home/shanghaixia/miniconda3/bin/fastqc
picard=/home/shanghaixia/bin/picard.jar
cutadapt=/home/shanghaixia/miniconda3/bin/cutadapt
trim=/home/shanghaixia/bin/Trimmomatic-0.39/trimmomatic-0.39.jar
adapter=/home/shanghaixia/bin/Trimmomatic-0.39/adapters
qualimap=/home/shanghaixia/bin/qualimap_v2.2.1/qualimap
### fastp去除接头和结尾
ID1=${1}-H3K27me3_L1_Q0181W0150

### mapping
cores=5

## post-alignment QC
## step 1: using picard tools for sorting and index aligned.bam file
${java} -jar ${picard} SortSam -I ${tar}/${ID}_bowtie2.sam -O ${tar}/${ID}_bowtie2.sorted.sam -SO coordinate --CREATE_INDEX TRUE


## step 2: using picard tool to mark duplicate
${java} -jar ${picard} MarkDuplicates I=${tar}/${ID}_bowtie2.sorted.sam O=${tar}/${ID}_bowtie2.sorted.sorted.dupMarked.sam METRICS_FILE=${tar}/${ID}_picard.dupMark.txt


## step 3: using picard tool to remove duplicate
${java} -jar ${picard} MarkDuplicates I=${tar}/${ID}_bowtie2.sorted.sam O=${tar}/${ID}_bowtie2.sorted.rmdup.sam REMOVE_DUPLICATES=true METRICS_FILE=${tar}/${ID}_picard.rmDup.txt

## step 4: filtering 

## Filter and keep the mapped read pairs
${samtools} view -bS -F 0x04 ${tar}/${ID}_bowtie2.sorted.rmdup.sam > ${tar}/${ID}_bowtie2.sorted.rmdup.bam

## step 4: remove  mitochondrial reads
${samtools} view -h ${tar}/${ID}_bowtie2.sorted.rmdup.bam | grep -v chrM | ${samtools} sort -O bam -o ${tar}/${ID}_bowtie2.sorted.rmdup.mapped.rmChrM.bam -T .

## step 5: smatools index
${samtools} index ${tar}/${ID}_bowtie2.sorted.rmdup.mapped.rmChrM.bam
${qualimap} bamqc -bam ${tar}/${ID}_bowtie2.sorted.rmdup.mapped.rmChrM.bam -gd hg19 -outdir . -outfile ${ID}.qualimap.report -outformat html
bamCoverage --bam ${tar}/${ID}_bowtie2.sorted.rmdup.mapped.rmChrM.bam -o ${tar}/${ID}_bowtie2.sorted.rmdup.mapped.rmChrM.Norm_new.bw --binSize 50 --normalizeUsing CPM --effectiveGenomeSize 2864785220 --extendReads


