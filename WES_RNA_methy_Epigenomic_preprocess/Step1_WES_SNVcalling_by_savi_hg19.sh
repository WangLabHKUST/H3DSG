#!/bin/bash
#SBATCH --mail-user=xxx #Update your email address
#SBATCH --mail-type=end
#SBATCH -p x-gpu-share
#SBATCH -N 1 -n 5


cd ${tar_dir}
export PATH="/home/zmoad/Software/miniconda2/install/bin:/scratch/PI/jgwang/zmoad/Software/savihg19/snpEff_v4_1c:$PATH" &&\
export LD_LIBRARY_PATH="/home/zmoad/Software/miniconda2/install/lib/"

PID=${1} ##PID
N=${PID}_B.sort.MD.bam ##blood sample
P=${PID}_T.sort.MD.bam ##tumor sample
nN='PID_blood' ##blood sample
nP='PID_tumor' ##tumor sample ID

dir=/scratch/PI/jgwang/haixia/projects/SP_WES/bam/hg19/new_hg19_map
python=/home/zmoad/Software/miniconda2/install/bin/python
S=/scratch/PI/jgwang/zmoad/Software/savihg19/SAVI-master/savi.py
V=/scratch/PI/jgwang/zmoad/Database/savihg19/vcf
R=/scratch/PI/jgwang/jtangbd/reference/hg19bwa/hg19.fa


mkdir /scratch/PI/jgwang/haixia/projects/SP_WES/savi_results/savi_hg19/${PID}

OUTPUT=/scratch/PI/jgwang/haixia/projects/SP_WES/savi_results/savi_hg19/${PID}
 

for i in {1..22} {X,Y}
do
	mkdir -p ${OUTPUT}/${i}
    ${S} --bams ${dir}/${N},${dir}/${P} \
		--names ${nN},${nP} \
		--ref ${R} \
		--region chr${i} -v \
		--outputdir ${OUTPUT}/${i} \
		--annvcf ${V}/cbio.fix.sort.vcf,${V}/meganormal186TCGA.fix.sort.vcf,${V}/CosmicCodingMuts.v72.May52015.jw.vcf,${V}/CosmicNonCodingVariants.v72.May52015vcf,${V}/All_20150605.vcf,${V}/219normals.cosmic.hitless100.noExactMut.mutless5000.all_samples.vcf 1> ${OUTPUT}/${i}/err.out 2>${OUTPUT}/${i}/err.log
done
wait

echo SAVI done!

echo Start merging...

O=${OUTPUT}/${nP}_vs_${nN}.report.coding.PDfilter.txt
head -n 1 ${OUTPUT}/1/report/report.coding.PDfilter.txt > ${O}

for i in {1..22} {X,Y} 
do
        F=${OUTPUT}/${i}/report/report.coding.PDfilter.txt
        if [ -e ${F} ]; then
                tail -n +2 ${F} >> ${O}
        else
                echo ${F} does not exist
        fi
done
