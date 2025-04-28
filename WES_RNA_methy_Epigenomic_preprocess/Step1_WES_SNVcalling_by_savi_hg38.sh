#!/bin/bash
#SBATCH --mail-user=shanghaixia@ust.hk #Update your email address
#SBATCH --mail-type=end
#SBATCH -p cpu-share
#SBATCH -N 1 -n 20

export PATH="/home/zmoad/Software/miniconda2/install/bin:/home/zmoad/Software/snpEff/snpEff:$PATH" &&\
export LD_LIBRARY_PATH="/home/zmoad/Software/miniconda2/install/lib/"
ID=P621820
savi_res_dir=/scratch/PI/jgwang/haixia/projects/SP_WES/savi_results/savi_hg38
cd ${savi_res_dir}

mkdir ${ID}
bam_dir=/scratch/PI/jgwang/haixia/projects/SP_WES/bam/hg38/hg38_new_map

Output_Dir=${savi_res_dir}/${ID}

if [ ! -d "${Output_Dir}" ]; then
mkdir ${Output_Dir}
fi &&\

V=/d3/scratch/PI/jgwang/zmoad/Database/SAVI_VCF &&\
S=/home/zmoad/Software/savi/SAVI3 &&\
#R=/home/zmoad/database/hg38chr/GRCh38.d1.vd1.fa &&\
R=/scratch/PI/jgwang/jtangbd/reference/TCGA_GRCh38_bwa_ref/GRCh38.d1.vd1.fa &&\



/home/zmoad/Software/miniconda2/install/bin/python /home/zmoad/Software/savi/SAVIgerm/SAVIgerm/savi_SONG_allregion_printCommand_allPDfilterSomatic_WithNullAnno.py --step 5 --conf 1e-3 --presence 1e-5 --noclean --bams ${bam_dir}/${ID}_B.sort.MD.bam,${bam_dir}/${ID}_T.sort.MD.bam --names ${ID}_B,${ID}_T --ann hg38_Song -v --superverbose --ref ${R} -v --outputdir ${Output_Dir} --annvcf ${V}/Song_hg38_Cosmic.sorted.vcf.gz,${V}/Song_hg38_ClinVar.sorted.vcf.gz,${V}/Song_hg38_PancancerGermline.sorted.vcf.gz,${V}/Song_hg38_cbio.vcf.gz,${V}/Song_hg38_All_SNP.sorted.vcf.gz,${V}/Song_hg38_TOPMED.sorted.vcf.gz,${V}/IGCTfamB.sorted.vcf.gz,${V}/IGCTpatB.sorted.vcf.gz,${V}/Song_hg38_GATKnormal.vcf.gz,${V}/Song_hg38_HGDP_1KG.sorted.vcf.gz,${V}/Song_hg38_MuTectPON.vcf.gz,${V}/Song_hg38_dbSNP_ALFA.sorted.vcf.gz,${V}/Song_hg38_gnomAD.sorted.vcf.gz,${V}/Song_hg38_219normals.vcf.gz,${V}/Song_China_map.sorted.vcf.gz,${V}/Song_hg38_ExAC.sorted.vcf.gz,${V}/Song_hg38_NGDC_SNP.vcf.gz,${V}/Song_hg38_meganormal.vcf.gz 1>>run.savi.hg38_1.sh.e 2>>run.savi.hg38_1.sh.o &&\
echo "jobs completed">>run.savi.hg38_1.sh.sign
