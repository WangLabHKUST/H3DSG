#!/bin/bash
#SBATCH --mail-user=shanghaixia@ust.hk #Update your email address
#SBATCH --mail-type=end
#SBATCH -p amd
#SBATCH -N 1 -n 20

export LD_LIBRARY_PATH="/home/zmoad/Software/miniconda2/install/lib/"

ID=P621820
tar_dir=/scratch/PI/jgwang/haixia/projects/SP_WES/cnv_results/cnv_hg38
cd ${tar_dir}
mkdir ${ID}

bam_dir=/scratch/PI/jgwang/haixia/projects/SP_WES/bam/hg38/hg38_new_map

cd ${ID} &&\

if [ ! -f "C3N_03188_Tumor-VS-C3N_03188_Normal.pileup" ];then
	# rm C3N_03188_Tumor-VS-C3N_03188_Normal.pileup
	echo neddpileup
	sh /home/zmoad/Software/facets/snp-pileup.hg38.WGS.sh ${bam_dir}/${ID}_B.sort.MD.bam ${bam_dir}/${ID}_T.sort.MD.bam ${ID}_T-VS-${ID}_B
fi &&\



#eval "$(conda shell.bash hook)" &&\
#conda activate R4.0.2  &&\

#Rscript /home/zmoad/Software/facets/run.facets.hg38.R ${ID}_T-VS-${ID}_B 25 &&\
/home/zmoad/Software/miniconda3/install/bin/Rscript /home/zmoad/Software/facets/run.facets.hg38.R ${ID}_T-VS-${ID}_B 50 &&\
/home/zmoad/Software/miniconda3/install/bin/Rscript /home/zmoad/Software/facets/run.facets.hg38.R ${ID}_T-VS-${ID}_B 100 &&\
/home/zmoad/Software/miniconda3/install/bin/Rscript /home/zmoad/Software/facets/run.facets.hg38.R ${ID}_T-VS-${ID}_B 150 &&\

# Rscript /home/zmoad/Software/facets/run.facets.hg38.R C3N_03188_Tumor-VS-C3N_03188_Normal 150 &&\
# Rscript /home/zmoad/Software/facets/run.facets.hg38.R C3N_03188_Tumor-VS-C3N_03188_Normal 200 &&\
# Rscript /home/zmoad/Software/facets/run.facets.hg38.R C3N_03188_Tumor-VS-C3N_03188_Normal 300 &&\
# Rscript /home/zmoad/Software/facets/run.facets.hg38.R C3N_03188_Tumor-VS-C3N_03188_Normal 400 &&\
# Rscript /home/zmoad/Software/facets/run.facets.hg38.R C3N_03188_Tumor-VS-C3N_03188_Normal 500 &&\
# Rscript /home/zmoad/Software/facets/run.facets.hg38.R C3N_03188_Tumor-VS-C3N_03188_Normal 1000 &&\
# Rscript /home/zmoad/Software/facets/run.facets.hg38.R C3N_03188_Tumor-VS-C3N_03188_Normal 1500 &&\
# Rscript /home/zmoad/Software/facets/run.facets.hg38.R C3N_03188_Tumor-VS-C3N_03188_Normal 2000 &&\
# Rscript /home/zmoad/Software/facets/run.facets.hg38.R C3N_03188_Tumor-VS-C3N_03188_Normal 5000 &&\
# Rscript /home/zmoad/Software/facets/run.facets.hg38.R C3N_03188_Tumor-VS-C3N_03188_Normal 10000 &&\
# Rscript /home/zmoad/Software/facets/run.facets.hg38.R C3N_03188_Tumor-VS-C3N_03188_Normal 50000 &&\
# Rscript /home/zmoad/Software/facets/run.facets.hg38.R C3N_03188_Tumor-VS-C3N_03188_Normal 100000 &&\


echo "jobs completed">>run.facets.hg38_1.sh.sign
