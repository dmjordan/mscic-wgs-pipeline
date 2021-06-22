# Daniel Jordan's adapter to coloc2
# to be run from Snakemake

library(readr)
library(stringr)
library(magrittr)
library(dplyr)
here::i_am("do_coloc2.R")
source(here::here("functions_coloc_likelihood_summary_integrated.R"))

eqtl_df <- read_tsv(snakemake@input[["eqtl"]])
split_id <- str_split_fixed(eqtl_df$variant_id, "_", 5)
eqtl_df %<>% transmute(SNPID=variant_id,
                       CHR=split_id[,1],
                       POS=as.integer(split_id[,2]),
                       A1=split_id[,4], A2=split_id[,3],
                       BETA=slope, SE=slope_se, F=maf, PVAL=pval_beta,
                       ProbeID=gene_id, N=948)

type <- if ("n.cases" %in% names(assoc_df)) "cc" else "quant"
assoc_df <- formatColoc(snakemake@input[["assoc"]], type=type, eqtl=FALSE)

coloc_results <- coloc.eqtl.biom(eqtl_df, assoc_df, outfolder="coloc2",
                                 prefix=snakemake@params[["prefix"]], cores=snakemake@resources[["cpus"]],
                                 save.coloc.output=TRUE, match_snpid=FALSE)

write_tsv(coloc_results, paste0("coloc2/", prefix, "_full_table.txt"))

NULL