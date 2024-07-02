#!/bin/bash
#SBATCH --mail-user=shanghaixia@ust.hk #Update your email address
#SBATCH --mail-type=end
#SBATCH -p himem-share
#SBATCH -N 1 -n 20

conda activate arriba_Bioconda
arriba=/home/shanghaixia/miniconda3/envs/arriba_Bioconda/bin/arriba
cd /scratch/PI/jgwang/haixia/projects/SP_RNA/gene_fusion_arriba_verion2_4/arriba_gene_fusion_results

dir=/scratch/PI/jgwang/haixia/projects/SP_RNA/fastq
R3=${1}
R1=${R3}_1.fq.gz
R2=${R3}_2.fq.gz

mkdir ${R3}_arriba_fusion
cd ${R3}_arriba_fusion

arriba_path=/home/shanghaixia/miniconda3/envs/arriba_Bioconda/var/lib/arriba

STAR_INDEX_DIR=/scratch/PI/jgwang/haixia/projects/SP_RNA/gene_fusion_arriba_verion2_4/reference/STAR_index_hs37d5viral_GENCODE19 #STAR_index_hs37d5_GENCODE19
ANNOTATION_GTF=${arriba_path}/GENCODE19.gtf 
ASSEMBLY_FA=${arriba_path}/hs37d5viral.fa 
BLACKLIST_TSV=${arriba_path}/blacklist_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz
KNOWN_FUSIONS_TSV=${arriba_path}/known_fusions_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz 
TAGS_TSV=${arriba_path}/known_fusions_hg19_hs37d5_GRCh37_v2.4.0.tsv.gz 
PROTEIN_DOMAINS_GFF3=${arriba_path}/protein_domains_hg19_hs37d5_GRCh37_v2.4.0.gff3 
THREADS=20
READ1=${dir}/${R1}
READ2=${dir}/${R2}

STAR \
    --runThreadN $THREADS \
    --genomeDir $STAR_INDEX_DIR --genomeLoad NoSharedMemory \
    --readFilesIn $READ1 $READ2 --readFilesCommand zcat \
    --outStd BAM_Unsorted --outSAMtype BAM Unsorted --outSAMunmapped Within --outBAMcompression 0 \
    --outFilterMultimapNmax 50 --peOverlapNbasesMin 10 --alignSplicedMateMapLminOverLmate 0.5 --alignSJstitchMismatchNmax 5 -1 5 5 \
    --chimSegmentMin 10 --chimOutType WithinBAM HardClip --chimJunctionOverhangMin 10 --chimScoreDropMax 30 \
    --chimScoreJunctionNonGTAG 0 --chimScoreSeparation 1 --chimSegmentReadGapMax 3 --chimMultimapNmax 50 | tee Aligned.out.bam |
$arriba \
    -x /dev/stdin \
   -o ${R3}_fusions.tsv -O ${R3}_fusions.discarded.tsv \
    -a $ASSEMBLY_FA -g $ANNOTATION_GTF \
    -b $BLACKLIST_TSV -k $KNOWN_FUSIONS_TSV