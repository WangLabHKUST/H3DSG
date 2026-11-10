
# Tool locations: set these for your system before running.
SOFTWARE_DIR=${SOFTWARE_DIR:-./software}
DATABASE_DIR=${DATABASE_DIR:-./database}
PROJECT_DIR=${PROJECT_DIR:-.}
REFERENCE_DIR=${REFERENCE_DIR:-./reference}

### 0220
./run_MutSigCV.sh ${PROJECT_DIR}/desk/MutSigCV_1.41/v901/ SP_Wild_21_split_new0220.maf exome_full192.coverage.txt gene.covariates.txt SP_Wild_21_split_new0220 mutation_type_dictionary_file.txt chr_files_hg19

