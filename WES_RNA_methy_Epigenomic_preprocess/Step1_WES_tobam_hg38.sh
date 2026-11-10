#!/bin/bash
#SBATCH --mail-user=your_email@example.org #Update your email address
#SBATCH --mail-type=end
#SBATCH -p amd
#SBATCH -N 1 -n 5

# Tool locations: set these for your system before running.
SOFTWARE_DIR=${SOFTWARE_DIR:-./software}
DATABASE_DIR=${DATABASE_DIR:-./database}
CONDA_DIR=${CONDA_DIR:-./miniconda3}
TOOL_BIN_DIR=${TOOL_BIN_DIR:-./bin}
TMP_DIR=${TMP_DIR:-./tmp}
PROJECT_DIR=${PROJECT_DIR:-.}
REFERENCE_DIR=${REFERENCE_DIR:-./reference}

R3=P621820_T
R1=${R3}_fp.R1.fq.gz 
R2=${R3}_fp.R2.fq.gz 


cd ${PROJECT_DIR}/SP_WES/bam/hg38/hg38_new_map
dir=${PROJECT_DIR}/SP_WES/fastq
out=${PROJECT_DIR}/SP_WES/bam/hg38/hg38_new_map
picard=${TOOL_BIN_DIR}/picard.jar
Ref1=${REFERENCE_DIR}/TCGA_GRCh38_bwa_ref/GRCh38.d1.vd1.fa
bwa=${CONDA_DIR}/bin/bwa
java=${CONDA_DIR}/bin/java
samtools=${CONDA_DIR}/bin/samtools

${bwa} mem -t 24 -M ${Ref1} ${dir}/${R1} ${dir}/${R2} | ${samtools} view -bSh -@ 8 - > ${out}/${R3}.bam
${samtools} sort -@ 24 ${out}/${R3}.bam ${out}/${R3}.sort
${java} -Djava.io.tmpdir=${TMP_DIR} -Xmx8g -jar ${picard} MarkDuplicates INPUT=${out}/${R3}.sort.bam OUTPUT=${out}/${R3}.sort.MD.bam METRICS_FILE=${out}/${R3}.sort.MD.bam.txt ASSUME_SORTED=true REMOVE_DUPLICATES=true VALIDATION_STRINGENCY=LENIENT
${samtools} index ${out}/${R3}.sort.MD.bam
rm ${out}/${R3}.bam ${out}/${R3}.sort.bam ${out}/${R3}.sort.MD.bam.txt
