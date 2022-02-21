library(GENESIS)
library(GWASTools)
library(tibble)
library(magrittr)
library(purrr)
library(stringr)
library(readr)
library(dplyr)

run_pcair <- function(gdsfile, kingfile, out_prefix) {
    rdsfile <- paste(out_prefix, "PCAir", "RDS", sep=".")
    txtfile <- paste(out_prefix, "PCAir", "txt", sep=".")

    # load genotype data
    geno <- GdsGenotypeReader(filename=gdsfile)
    genoData <- GenotypeData(geno)

    # load KING result
    KINGmat <- kingToMatrix(kingfile, estimator="Kinship", sample.include=getVariable(genoData, "sample.id"))

    # run PC-AIR
    pcair_result <- pcair(genoData, kinobj=KINGmat, divobj=KINGmat)

    # Plot first 10 PCs
    pc_table <- as_tibble(pcair_result$vectors)
    pc_table %>% rename_all(~ str_replace(.x, 'V', 'PC')) %>%
                 add_column(Subject_ID=pcair_result$sample.id, .before=1) -> pc_table
    saveRDS(pcair_result, file=rdsfile)
    pc_table %>% select(c(Subject_ID, starts_with("PC"))) %>% write_delim(path=txtfile)
    pc_table
}

run_pcair(snakemake@input[["gds"]], snakemake@input[["king"]], snakemake@params[["output_stem"]])

NULL
