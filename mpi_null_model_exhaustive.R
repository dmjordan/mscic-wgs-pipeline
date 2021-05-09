suppressPackageStartupMessages(library(doMPI))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GWASTools))
suppressPackageStartupMessages(library(GENESIS))

cl <- startMPIcluster()
registerDoMPI(cl)

args <- commandArgs(trailingOnly=TRUE)
file_prefix <- args[[1]]
endpoint <- args[[2]]

read_csv("/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/jordad05/625_Samples.cohort.QC_filtered.sample_matched.age_flowcell_PCAir_dmatrix.csv") %>%
  rename(scanID=X1) %>% mutate(race_factor = factor(race_factor)) -> clinical_table
scan_annot <- ScanAnnotationDataFrame(as.data.frame(clinical_table))  # somehow GWASTools doesn't recognize tibble columns?

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

map(1:10, ~paste0("pc", 1:.x)) %>% prepend(list(character())) -> pcs
map(1:length(non_pc_covars), ~ combn(non_pc_covars, .x, simplify=FALSE)) %>%
  flatten %>% cross2(pcs) %>% map(unlist) -> all_covar_combinations

family <- if (clinical_table %>% filter(.[[endpoint]] != 0 & .[[endpoint]] != 1) %>% tally > 0) "gaussian" else "binomial"

pcrel_result <- readRDS(paste(file_prefix, "PCRelate", "RDS", sep="."))
grm <- pcrelateToMatrix(pcrel_result)

ptm <- proc.time()
models <- foreach (covars=all_covar_combinations,
         .packages="GENESIS",
         .errorhandling="remove") %dopar% {
        model <- fitNullModel(scan_annot, outcome=endpoint,
                              covars=covars,
                              cov.mat=grm, family=family)
        model$covars <- covars
        model
    }
exec_time <- proc.time() - ptm
cat(length(models), "of", length(all_covar_combinations), "covar combinations tested in", exec_time[["elapsed"]], "seconds\n")
map_dbl(models, "AIC") %>% compact %>% which.min %>% pluck(models, .) -> model
cat("Best combination:", model$covars, "\n")

saveRDS(model, paste(file_prefix, endpoint, "null", "RDS", sep="."))