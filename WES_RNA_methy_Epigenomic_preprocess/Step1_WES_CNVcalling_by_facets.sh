#!/bin/bash
#SBATCH --mail-user=XXX #Update your email address
#SBATCH --mail-type=end
#SBATCH -p x-gpu
#SBATCH -N 1 -n 5


# Tool locations: set these for your system before running.
SOFTWARE_DIR=${SOFTWARE_DIR:-./software}
DATABASE_DIR=${DATABASE_DIR:-./database}
CONDA_DIR=${CONDA_DIR:-./miniconda3}
TOOL_BIN_DIR=${TOOL_BIN_DIR:-./bin}
TMP_DIR=${TMP_DIR:-./tmp}
PROJECT_DIR=${PROJECT_DIR:-.}
REFERENCE_DIR=${REFERENCE_DIR:-./reference}

perl ${SOFTWARE_DIR}/facets19/pileup/snp-pileup.WithChr.pl ${PROJECT_DIR}/SP_WES/bam/hg19/new_hg19_map/chenyulin-QX.sort.MD.bam ${PROJECT_DIR}/SP_WES/bam/hg19/new_hg19_map/chenyulin_L1.sort.MD.bam chenyulin_L1-VS-chenyulin_QX.mpileup &&\

export LD_LIBRARY_PATH="${SOFTWARE_DIR}/miniconda2/install/lib/"

${SOFTWARE_DIR}/miniconda3/install/bin/Rscript ${SOFTWARE_DIR}/facets19/demo/run.facets.hg19.R chenyulin_L1-VS-chenyulin_QX 150
${SOFTWARE_DIR}/miniconda3/install/bin/Rscript ${SOFTWARE_DIR}/facets19/demo/run.facets.hg19.R chenyulin_L1-VS-chenyulin_QX 200



