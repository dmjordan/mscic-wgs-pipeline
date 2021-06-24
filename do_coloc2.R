# Daniel Jordan's adapter to coloc2
# to be run from Snakemake

library(readr)
library(stringr)
library(magrittr)
library(dplyr)
snakemake@source("functions_coloc_likelihood_summary_integrated.R")
snakemake@source("optim_function.R")

cat("loading eqtls from", snakemake@input[["eqtl"]], "\n")
eqtl_df <- read_tsv(snakemake@input[["eqtl"]])
cat("loaded", nrow(eqtl_df), "lines\n")
split_id <- str_split_fixed(eqtl_df$variant_id, "_", 5)
eqtl_df %<>% transmute(SNPID=variant_id,
                       CHR=split_id[,1],
                       POS=as.integer(split_id[,2]),
                       A1=split_id[,4], A2=split_id[,3],
                       BETA=slope, SE=slope_se, F=maf, PVAL=pval_beta,
                       ProbeID=gene_id, N=948)
cat(nrow(eqtl_df), "lines remain after transform\n")

cat("loading association results from from", snakemake@input[["assoc"]], "\n")
con <- file(snakemake@input[["assoc"]],"r")
first_line <- readLines(con,n=1)
close(con)
type <- if (str_detect(first_line, "n\\.cases")) "cc" else "quant"
cat("detected a", type, "phenotype\n")
assoc_df <- formatColoc(snakemake@input[["assoc"]], type=type, eqtl=FALSE)
cat("loaded", nrow(assoc_df), "lines\n")
assoc_df$type <- type
coloc_results <- coloc.eqtl.biom(eqtl_df, assoc_df, outfolder="coloc2",
                                 prefix=snakemake@params[["prefix"]], cores=snakemake@resources[["cpus"]],
                                 match_snpid=FALSE)

write_tsv(coloc_results, paste0("coloc2/", snakemake@params[["prefix"]], ".full_table.txt"))

NULL
