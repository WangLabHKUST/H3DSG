#!/bin/bash
#SBATCH --mail-user=shanghaixia@ust.hk #Update your email address
#SBATCH --mail-type=end
#SBATCH -p cpu-share
#SBATCH -N 1 -n 20
dir=/scratch/PI/jgwang/haixia/projects/SP_ATAC_seq/fq
ID=${1}
cd /scratch/PI/jgwang/haixia/projects/SP_ATAC_seq_new_mapping
mkdir ${ID}_trim
cd ${ID}_trim
tar=/scratch/PI/jgwang/haixia/projects/SP_ATAC_seq_new_mapping/${ID}_trim
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
trim=/home/shanghaixia/bin/Trimmomatic-0.39/trimmomatic-0.39.jar
adapter=/home/shanghaixia/bin/Trimmomatic-0.39/adapters
### trim 去除接头
${java} -jar ${trim} PE ${dir}/${ID}_R1.fq.gz ${dir}/${ID}_R2.fq.gz ${dir}/${ID}_R1_Trim_paired.fq.gz ${dir}/${ID}_R1_Trim_unpaired.fq.gz ${dir}/${ID}_R2_Trim_paired.fq.gz ${dir}/${ID}_R2_Trim_unpaired.fq.gz ILLUMINACLIP:${adapter}/NexteraPE-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads 10

### fastqc
${fastqc} -t 2 -o ${dir}/qc/ -f fastq ${dir}/${ID}_R1_Trim_paired.fq.gz
${fastqc} -t 2 -o ${dir}/qc/ -f fastq ${dir}/${ID}_R2_Trim_paired.fq.gz
### mapping

${bowtie2} --threads 8 -x ${Ref1} -q -1 ${dir}/${ID}_R1_Trim_paired.fq.gz -2 ${dir}/${ID}_R2_Trim_paired.fq.gz -S ${tar}/${ID}_new_trim.sam
${samtools} view -@ 8 -bS ${tar}/${ID}_new_trim.sam > ${tar}/${ID}_new_trim.bam
## sort bam file
${samtools} sort -@ 8 ${tar}/${ID}_new_trim.bam ${tar}/${ID}_new_trim.sorted.bam
## index bam file
${samtools} index ${tar}/${ID}_new_trim.sorted.bam.bam
## clean sam file 
rm ${tar}/${ID}_new_trim.sam
rm ${tar}/${ID}_new_trim.bam.bam
## identify and remove reads aligning to Mt and Chloroplast

## call the BAM file without chloroplastic and mitochondrial alignments as
${samtools} idxstats ${tar}/${ID}_new_trim.sorted.bam.bam | cut -f1 | grep -v Mt | grep -v Pt | xargs samtools view -@ 7 -b ${tar}/${ID}_new_trim.sorted.bam.bam > ${tar}/${ID}_new_trim.sorted.noorg.bam
${samtools} index ${tar}/${ID}_new_trim.sorted.noorg.bam


rm ${tar}/${ID}_new_trim.sorted.bam.bam
rm ${tar}/${ID}_new_trim.sorted.bam.bam.bai
