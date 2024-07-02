#!/bin/bash
#SBATCH --mail-user=XXX #Update your email address
#SBATCH --mail-type=end
#SBATCH -p x-gpu
#SBATCH -N 1 -n 5


perl /scratch/PI/jgwang/zmoad/Software/facets19/pileup/snp-pileup.WithChr.pl /scratch/PI/jgwang/haixia/projects/SP_WES/bam/hg19/new_hg19_map/chenyulin-QX.sort.MD.bam /scratch/PI/jgwang/haixia/projects/SP_WES/bam/hg19/new_hg19_map/chenyulin_L1.sort.MD.bam chenyulin_L1-VS-chenyulin_QX.mpileup &&\

export LD_LIBRARY_PATH="/home/zmoad/Software/miniconda2/install/lib/"

/home/zmoad/Software/miniconda3/install/bin/Rscript /scratch/PI/jgwang/zmoad/Software/facets19/demo/run.facets.hg19.R chenyulin_L1-VS-chenyulin_QX 150
/home/zmoad/Software/miniconda3/install/bin/Rscript /scratch/PI/jgwang/zmoad/Software/facets19/demo/run.facets.hg19.R chenyulin_L1-VS-chenyulin_QX 200



