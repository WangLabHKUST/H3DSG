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

cd ${PROJECT_DIR}/SP_RNA/bamhg38
dir=${PROJECT_DIR}/SP_RNA/batch1_fastq
out=${PROJECT_DIR}/SP_RNA/bamhg38

cat ${dir}/batch1.fq.STAR.txt | while read line
do
        echo ${line}
        RR=$line
        arr=($RR)
        R1=${arr[0]}
        R2=${arr[1]}
        R3=${arr[2]}
        mkdir -p ${R3}.star
        ulimit -n 4096
        STAR    --runThreadN 20   --genomeDir ${PROJECT_DIR}/SP_RNA/hg38_genome_gtf/STAR_index_hg38_genecode38 --outFileNamePrefix ${R3}.star/${R3} \
        --readFilesIn ${dir}/${R1} ${dir}/${R2} \
        --readFilesCommand zcat      --limitBAMsortRAM 0   \
        --outSAMtype BAM   SortedByCoordinate      --outSAMstrandField intronMotif   \
        --outSAMattributes NH   HI   NM   MD   AS   XS      --outSAMunmapped Within  \
        --outSAMheaderHD @HD   VN:1.4      \
        --outFilterMultimapNmax 20   --outFilterMultimapScoreRange 1   \
        --outFilterScoreMinOverLread 0.33   --outFilterMatchNminOverLread 0.33   \
        --outFilterMismatchNmax 10   --alignIntronMax 500000   \
        --alignMatesGapMax 1000000   --alignSJDBoverhangMin 1   --sjdbOverhang 100   --sjdbScore 2
        samtools index ${R3}.star/${R3}.Aligned.sortedByCoord.out.bam
done

