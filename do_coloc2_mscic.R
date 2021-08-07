
# Daniel Jordan's adapter to coloc2
# to be run from Snakemake

library(readr)
library(stringr)
library(magrittr)
library(dplyr)
snakemake@source("functions_coloc_likelihood_summary_integrated.R")
snakemake@source("optim_function.R")

cat("loading eqtls from", snakemake@input[["eqtl"]], "\n")
eqtl_df <- read_delim(snakemake@input[['eqtl']], delim=" ",
                    col_names=c("ProbeID", "probeChr", "probeStart", "probeEnd", "strand", "NVariants", "distToTopVar", "SNPID", "CHR", "POS", "varEnd", "PVAL", "BETA", "flag"),
                    col_types=list(ProbeID=col_character(), probeChr=col_character(), strand=col_character(), SNPID=col_character(), CHR=col_character(), .default=col_double()))
cat("loaded", nrow(eqtl_df), "lines\n")

cat("loading allele info from", snakemake@input[['bim']], '\n')
bimfile <- read_tsv(snakemake@input[['bim']],
                    col_names=c("CHR", "SNPID", "CM", "POS", "A1", "A2"))

eqtl_df %>% left_join(bimfile, by="SNPID", suffix=c("", "_y")) %>%
  select(SNPID, CHR, POS, A1, A2, ProbeID, BETA, PVAL) -> eqtl_df
cat(nrow(eqtl_df), "rows remain after join\n")

cat("reading individuals list from", snakemake@input[["fam"]], "\n")
read_delim(snakemake@input[['fam']], delim=' ', col_names=FALSE) %>%
  count %>% pluck('n') -> N
cat("counted", N, 'individuals\n')
eqtl_df$N <- N

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
                                 prefix=snakemake@params[["prefix"]],  # cores=snakemake@resources[["cpus"]],
                                 useBETA=FALSE, save.coloc.output=TRUE, match_snpid=FALSE)

if (!file.exists("coloc2")) dir.create("coloc2")
write_tsv(coloc_results, paste0("coloc2/", snakemake@params[["prefix"]], ".full_table.txt"))

NULL

