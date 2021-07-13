snakemake@source("rnaseq-gofigure.R")
#library(tidyverse)
select = dplyr::select

coloc2_results <- read_tsv(snakemake@input[[1]])
geneSig <- (coloc2_results$min.pval.biom < 0.05) & (coloc2_results$PP.H4.abf > 0.1)
names(geneSig) <- str_split_fixed(coloc2_results$ProbeID, fixed("."), n=2)[,1]


res_table <- run_topGO_Fisher(geneSig)

sig_res_table <- res_table %>%
    mutate(adj.P.Val = pmax(adj.P.Val, 10^-300, na.rm = TRUE)) %>%
    filter(adj.P.Val <= 0.05)
## Put in "standard input" format:
## https://gitlab.com/evogenlab/GO-Figure#standard-input
standard_input_table <- sig_res_table %>%
    select(GO.ID, adj.P.Val)
write_tsv(standard_input_table, snakemake@output[[1]], col_names = FALSE)
