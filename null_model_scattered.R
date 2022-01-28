suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GWASTools))
suppressPackageStartupMessages(library(GENESIS))

race <- c("WHITE", "BLACK", "HISPANIC", "ASIAN")

load_clinical_table <- function(dmatrix_path) {
  read_csv(dmatrix_path) %>%
    rename(scanID=Subject_ID) %>% mutate(race_factor = factor(race_factor))
}

get_covars <- function (clinical_table, index) {

  clinical_table %>% names %>% str_subset("^flowcell") %>% list %>%
    append(list("age",
                "age_squared",
                "sex",
                "age_sex",
                "recruitment_date",
                "is_hispanic",
                "race_factor",
                "multi_batch",
                "any_comorbidity")) -> non_pc_covars

  map(1:10, ~paste0(race, "_pc", 1:.x)) %>% prepend(list(character())) -> pcs
  map(1:length(non_pc_covars), ~ combn(non_pc_covars, .x, simplify=FALSE)) %>%
    flatten %>% cross2(pcs) %>% map(unlist) -> all_covar_combinations

  all_covar_combinations[[index]]
}

generate_null_model <- function(endpoint, scan_annot, pcrel_path, covars, sample.id) {
  family <- if (clinical_table %>% filter(.[[endpoint]] != 0 & .[[endpoint]] != 1) %>% tally > 0) "gaussian" else "binomial"

  pcrel_result <- readRDS(pcrel_path)
  grm <- pcrelateToMatrix(pcrel_result)
  tryCatch({
      model <- fitNullModel(scan_annot, outcome=endpoint,
                            covars=covars,
                            sample.id=sample.id, cov.mat=grm, family=family)
      model$covars <- covars
      model
  }, error=function(cond) list(error=cond$message, converged=FALSE))
}

clinical_table <- load_clinical_table(snakemake@input[[1]])
scan_annot <- ScanAnnotationDataFrame(as.data.frame(clinical_table))  # somehow GWASTools doesn't recognize tibble columns?
covars <- get_covars(clinical_table, snakemake@wildcards[["index"]])
sample.id <- tryCatch(read_tsv(snakemake@input[["indiv_list"]])[[1]],
                      error=function(cond) { getScanID(scan_annot) })
nullmod <- generate_null_model(snakemake@wildcards[["phenotype_untagged"]], scan_annot, snakemake@input[[2]], covars, sample.id)

saveRDS(nullmod, snakemake@output[[1]])

NULL
