#!/bin/bash
#SBATCH --mail-user=shanghaixia@ust.hk #Update your email address
#SBATCH --mail-type=end
#SBATCH -p amd
#SBATCH -N 1 -n 5

R3=P621820_T
R1=${R3}_fp.R1.fq.gz 
R2=${R3}_fp.R2.fq.gz 


cd /scratch/PI/jgwang/haixia/projects/SP_WES/bam/hg38/hg38_new_map
dir=/scratch/PI/jgwang/haixia/projects/SP_WES/fastq
out=/scratch/PI/jgwang/haixia/projects/SP_WES/bam/hg38/hg38_new_map
picard=/home/shanghaixia/bin/picard.jar
Ref1=/scratch/PI/jgwang/jtangbd/reference/TCGA_GRCh38_bwa_ref/GRCh38.d1.vd1.fa
bwa=/home/shanghaixia/miniconda3/bin/bwa
java=/home/shanghaixia/miniconda3/bin/java
samtools=/home/shanghaixia/miniconda3/bin/samtools

${bwa} mem -t 24 -M ${Ref1} ${dir}/${R1} ${dir}/${R2} | ${samtools} view -bSh -@ 8 - > ${out}/${R3}.bam
${samtools} sort -@ 24 ${out}/${R3}.bam ${out}/${R3}.sort
${java} -Djava.io.tmpdir=/home/shanghaixia/tmp -Xmx8g -jar ${picard} MarkDuplicates INPUT=${out}/${R3}.sort.bam OUTPUT=${out}/${R3}.sort.MD.bam METRICS_FILE=${out}/${R3}.sort.MD.bam.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
${samtools} index ${out}/${R3}.sort.MD.bam
rm ${out}/${R3}.bam ${out}/${R3}.sort.bam ${out}/${R3}.sort.MD.bam.txt
