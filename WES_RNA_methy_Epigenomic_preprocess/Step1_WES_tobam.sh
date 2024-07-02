#!/bin/bash
#SBATCH --mail-user=shanghaixia@ust.hk #Update your email address
#SBATCH --mail-type=end
#SBATCH -p x-gpu-share
#SBATCH -N 1 -n 5

R1=${1}_R1.fastq.gz 
R2=${1}_R2.fastq.gz
R3=${1}

cd /scratch/PI/jgwang/haixia/projects/SP_WES/bam/hg19/new_hg19_map
dir=/scratch/PI/jgwang/haixia/projects/SP_WES/fastq
out=/scratch/PI/jgwang/haixia/projects/SP_WES/bam/hg19/new_hg19_map
picard=/home/shanghaixia/bin/picard.jar
Ref1=/scratch/PI/jgwang/jtangbd/reference/hg19bwa/hg19.fa
bwa=/home/shanghaixia/miniconda3/bin/bwa
java=/home/shanghaixia/miniconda3/bin/java
samtools=/home/shanghaixia/miniconda3/bin/samtools

${bwa} mem -t 24 -M ${Ref1} ${dir}/${R1} ${dir}/${R2} | ${samtools} view -bSh -@ 8 - > ${out}/${R3}.bam
${samtools} sort -@ 24 ${out}/${R3}.bam ${out}/${R3}.sort
${java} -Djava.io.tmpdir=/home/shanghaixia/tmp -Xmx8g -jar ${picard} MarkDuplicates INPUT=${out}/${R3}.sort.bam OUTPUT=${out}/${R3}.sort.MD.bam METRICS_FILE=${out}/${R3}.sort.MD.bam.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
${samtools} index ${out}/${R3}.sort.MD.bam
rm ${out}/${R3}.bam ${out}/${R3}.sort.bam ${out}/${R3}.sort.MD.bam.txt
