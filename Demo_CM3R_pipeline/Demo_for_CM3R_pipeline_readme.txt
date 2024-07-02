1. Multi-omics data preparation
(a) gene expression: preprocessed by TPM, log2(), quantile normalization, min-max normlaization;
copy number variant: in the .seg mean format; CpG probes: beta values
(b) gene names with differential expression in the cohort
(c) CpG names with differential methylation CpG probes
2. Constrined_Ridge_regression (in Python environment)
Open 2_constrained_Ridge_regression/CNV1_cpg_3_TF_3_RegNetwork_Demo_E2F1_MYC.py in python
change the work directory
Run the CNV1_cpg_3_TF_3_RegNetwork_Demo_E2F1_MYC.py

3. Regulator identification 
Run 3_regualtor_identification/Regulator_identification_E2F1_MYC_Demo.R (in R environment)
the regulation network will be output in .csv
the master regulator will be output in the R code