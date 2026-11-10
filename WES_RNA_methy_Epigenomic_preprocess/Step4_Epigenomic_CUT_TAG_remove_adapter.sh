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

dir=${PROJECT_DIR}/SP_CUT_Tag
ID=${1}
cd ${PROJECT_DIR}/SP_CUT_TAG_new_mapping
#mkdir ${ID}_cut_tag_remove_adapter
cd ${ID}_cut_tag_remove_adapter
tar=${PROJECT_DIR}/SP_CUT_TAG_new_mapping/${ID}_cut_tag_remove_adapter
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
cutadapt=${CONDA_DIR}/bin/cutadapt
trim=${TOOL_BIN_DIR}/Trimmomatic-0.39/trimmomatic-0.39.jar
adapter=${TOOL_BIN_DIR}/Trimmomatic-0.39/adapters
qualimap=${TOOL_BIN_DIR}/qualimap_v2.2.1/qualimap
### fastp: remove adapters and trim read ends
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


