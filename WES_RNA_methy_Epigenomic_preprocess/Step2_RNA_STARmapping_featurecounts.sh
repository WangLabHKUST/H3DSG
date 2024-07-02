### star step 1: generating genome index files
srun -p x-gpu-share -N 1 -n 1 --pty bash
cd /scratch/PI/jgwang/haixia/projects/SP_RNA/ref
STAR --runMode genomeGenerate --genomeDir ref_star_hg19/ --genomeFastaFiles hg19.fa \
--sjdbOverhang 100 \
--sjdbGTFfile Homo_sapiens.GRCh37.75.chr.gtf --runThreadN 16


### star step 2: mapping reads to the genome
srun -p x-gpu -N 1 -n 1 --pty bash
cd /scratch/PI/jgwang/haixia/projects/SP_RNA/bam
dir=/scratch/PI/jgwang/haixia/projects/SP_RNA/fastq
out=/scratch/PI/jgwang/haixia/projects/SP_RNA/bam


cat ${dir}/SP_RNA_Batch_STAR.txt | while read line
do
        echo ${line}
        RR=$line
        arr=($RR)
        R1=${arr[0]}
        R2=${arr[1]}
        R3=${arr[2]}
        mkdir -p ${R3}.star
        ulimit -n 4096
        STAR    --runThreadN 20   --genomeDir /scratch/PI/jgwang/haixia/projects/SP_RNA/ref/ref_star_hg19 --outFileNamePrefix ${R3}.star/${R3} \
        --readFilesIn ${dir}/${R1} ${dir}/${R2} \
        --readFilesCommand zcat      --limitBAMsortRAM 0   \
        --outSAMtype BAM   SortedByCoordinate      --outSAMstrandField intronMotif   \
        --outSAMattributes NH   HI   NM   MD   AS   XS      --outSAMunmapped Within  \
        --outSAMheaderHD @HD   VN:1.4      \
        --outFilterMultimapNmax 20   --outFilterMultimapScoreRange 1   \
        --outFilterScoreMinOverLread 0.33   --outFilterMatchNminOverLread 0.33   \
        --outFilterMismatchNmax 10   --alignIntronMax 500000   \
        --alignMatesGapMax 1000000   --alignSJDBoverhangMin 1   --sjdbOverhang 50   --sjdbScore 2
        samtools index ${R3}.star/${R3}Aligned.sortedByCoord.out.bam
done


#### step 3: generate count table
srun -p cpu-share -N 1 -n 1 --pty bash
cd /scratch/PI/jgwang/haixia/projects/SP_RNA/bam

pid=${1}
I_bam=/scratch/PI/jgwang/haixia/projects/SP_RNA/bam/${pid}.star/${pid}Aligned.sortedByCoord.out.bam
O_out=/scratch/PI/jgwang/haixia/projects/SP_RNA/bam/${pid}.star/${pid}count
featureCounts -p -T 12 -B -t exon -g gene_name -a /scratch/PI/jgwang/haixia/projects/SP_RNA/ref/Homo_sapiens.GRCh37.75.chr.gtf -o $O_out"_I_hg19.txt" $I_bam



