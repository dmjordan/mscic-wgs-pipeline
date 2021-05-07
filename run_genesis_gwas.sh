#!/bin/bash
set -eo pipefail
ml plink

bsub -K < ../../scripts/WGS/lsf_submit_prep_files_for_pcair.sh
plink --bfile 625_Samples.cohort.QC_filtered.sample_matched.maf_and_ld --make-bed --out 625_Samples.cohort.QC_filtered.sample_matched.maf_and_ld.plink_recoded
conda run -n r-gwas Rscript --vanilla ../../scripts/WGS/genesis_pca.R
