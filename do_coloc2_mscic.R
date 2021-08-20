
# Daniel Jordan's adapter to coloc2

suppressPackageStartupMessages({
#library(doMPI)
library(foreach)
library(readr)
library(stringr)
library(magrittr)
library(tidyr)
library(dplyr)
library(purrr)
})

#cl <- startMPIcluster()
#registerDoMPI(cl)

#args <- commandArgs(trailingOnly=TRUE)
#script_dir <- args[[1]]
#eqtl_file <- args[[2]]
#sample_file <- args[[3]]
#assoc_file <- args[[4]]
#output_prefix <- args[[5]]

eqtl_file <- snakemake@input[["eqtl"]]
sample_file <- snakemake@input[["bim"]]
assoc_file <- snakemake@input[["assoc"]]
output_prefix <- snakemake@params[["prefix"]]

snakemake@source("functions_coloc_likelihood_summary_integrated.R")
snakemake@source("optim_function.R")

cat("loading eqtls from", eqtl_file, "\n")
eqtl_df <- read_delim(eqtl_file, delim=" ",
                    col_types=list(ProbeID=col_character(), SNPID=col_character(),
                                   CHR=col_character(), POS=col_integer(), A1=col_character(), A2=col_character(),
                                   BETA=col_number(), PVAL=col_number()))
cat("loaded", nrow(eqtl_df), "lines\n")

eqtl_df %<>% arrange(CHR, POS, A1, A2)

cat("reading individuals list from", sample_file, "\n")
read_tsv(sample_file, col_names=FALSE) %>%
  count %>% pluck('n') -> N
cat("counted", N, 'individuals\n')
eqtl_df$N <- N

cat("loading association results from ", sample_file, "\n")
con <- file(sample_file,"r")
first_line <- readLines(con,n=1)
close(con)
type <- if (str_detect(first_line, "n\\.cases")) "cc" else "quant"
cat("detected a", type, "phenotype\n")
assoc_df <- formatColoc(assoc_file, type=type, eqtl=FALSE)
cat("loaded", nrow(assoc_df), "lines\n")
assoc_df$type <- type
if (!file.exists("coloc2")) dir.create("coloc2")
#chunks <- eqtl_df %>% group_by(chunk=(row_number()-1) %/% (n()/length(cl))) %>%
                    nest %>% pull(data)
#coloc_results <- foreach(eqtl_df=chunks,
#        .combine=bind_rows) %dopar% {
    coloc_results <- coloc.eqtl.biom(eqtl_df, assoc_df, outfolder="coloc2",
                                     prefix=output_prefix,  # cores=snakemake@resources[["cpus"]],
                                     useBETA=FALSE, save.coloc.output=FALSE, match_snpid=FALSE)
#}

write_tsv(coloc_results, paste0("coloc2/", output_prefix, ".full_table.txt"))

#closeCluster(cl)
#mpi.quit()

NULL

