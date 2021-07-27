snakemake@source("rnaseq-gofigure.R")
#library(tidyverse)
select = dplyr::select

results <- read_delim(snakemake@input[[1]],
                      delim = if (snakemake@params[["data_type"]] == "spredixcan") "," else "\t")
geneSig <- if (snakemake@params[["data_type"]] == "coloc2") {
  (results$min.pval.biom <= 0.05) & (results$PP.H4.abf >= 0.1)
} else {
  results$pvalue <= 0.01
}
ensemblNames <- if (snakemake@params[["data_type"]] == "coloc2") results$ProbeID else results$gene
names(geneSig) <- str_split_fixed(ensemblNames, fixed("."), n=2)[,1]

res_table <- run_topGO_Fisher(geneSig)

sig_res_table <- res_table %>%
    mutate(adj.P.Val = pmax(adj.P.Val, 10^-300, na.rm = TRUE)) %>%
    filter(adj.P.Val <= 0.05)
## Put in "standard input" format:
## https://gitlab.com/evogenlab/GO-Figure#standard-input
standard_input_table <- sig_res_table %>%
    select(GO.ID, adj.P.Val)
write_tsv(standard_input_table, snakemake@output[[1]], col_names = FALSE)
