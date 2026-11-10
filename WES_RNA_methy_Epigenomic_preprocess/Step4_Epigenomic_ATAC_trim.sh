#!/bin/bash
#SBATCH --mail-user=your_email@example.org #Update your email address
#SBATCH --mail-type=end
#SBATCH -p cpu-share
#SBATCH -N 1 -n 20
# Tool locations: set these for your system before running.
SOFTWARE_DIR=${SOFTWARE_DIR:-./software}
DATABASE_DIR=${DATABASE_DIR:-./database}
CONDA_DIR=${CONDA_DIR:-./miniconda3}
TOOL_BIN_DIR=${TOOL_BIN_DIR:-./bin}
TMP_DIR=${TMP_DIR:-./tmp}
PROJECT_DIR=${PROJECT_DIR:-.}
REFERENCE_DIR=${REFERENCE_DIR:-./reference}

dir=${PROJECT_DIR}/SP_ATAC_seq/fq
ID=${1}
cd ${PROJECT_DIR}/SP_ATAC_seq_new_mapping
mkdir ${ID}_trim
cd ${ID}_trim
tar=${PROJECT_DIR}/SP_ATAC_seq_new_mapping/${ID}_trim
Ref1=${PROJECT_DIR}/SP_Chipseq/bowtie2_indexes/hg19

java=${CONDA_DIR}/bin/java
samtools=${CONDA_DIR}/bin/samtools
bowtie2=${CONDA_DIR}/bin/bowtie2
macs2=${CONDA_DIR}/bin/macs2
bioawk=${CONDA_DIR}/bin/bioawk
bedtools=${CONDA_DIR}/bin/bedtools
sort=/usr/bin/sort
bedClip=${CONDA_DIR}/bin/bedClip
bedGraphToBigWig=${CONDA_DIR}/bin/bedGraphToBigWig
fastp=${CONDA_DIR}/bin/fastp
fastqc=${CONDA_DIR}/bin/fastqc
picard=${TOOL_BIN_DIR}/picard.jar
trim=${TOOL_BIN_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar
adapter=${TOOL_BIN_DIR}/Trimmomatic-0.39/adapters
### trim: remove adapters
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
