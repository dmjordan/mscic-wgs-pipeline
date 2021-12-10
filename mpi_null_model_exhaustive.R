suppressPackageStartupMessages(library(doMPI))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(GWASTools))
suppressPackageStartupMessages(library(GENESIS))

cl <- startMPIcluster()
registerDoMPI(cl)

args <- commandArgs(trailingOnly=TRUE)
file_prefix <- args[[1]]
endpoint <- args[[2]]
subset_tag <- if (length(args) > 2) args[[3]] else ""
exclude_samples <- if (length(args) > 3) as.character(args[4:length(args)]) else character(0)

cat("parsed args:\n")
cat("file prefix:", file_prefix, fill=TRUE)
cat("endpoint:", endpoint, fill=TRUE)
cat("subset:", if (subset_tag == "") "<none>" else subset_tag, fill=TRUE)
cat("excluded samples:", if (length(exclude_samples) > 0) exclude_samples else "<none>", fill=TRUE)

read_csv("/sc/arion/projects/mscic1/data/covariates/clinical_data_deidentified_allsamples/jordad05/625_Samples.cohort.QC_filtered.sample_matched.age_flowcell_PCAir_dmatrix.csv") %>%
  rename(scanID=Subject_ID) %>% mutate(race_factor = factor(race_factor)) -> clinical_table
scan_annot <- ScanAnnotationDataFrame(as.data.frame(clinical_table))  # somehow GWASTools doesn't recognize tibble columns?
sample.id <- setdiff(getScanID(scan_annot), exclude_samples)

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

pcs <- paste0("pc", 1:10)
map(1:length(non_pc_covars), ~ combn(non_pc_covars, .x, simplify=FALSE)) %>%
  flatten %>% map(unlist) %>% map(~ c(.x, pcs)) -> all_covar_combinations

family <- if (clinical_table %>% filter(.[[endpoint]] != 0 & .[[endpoint]] != 1) %>% tally > 0) "gaussian" else "binomial"

pcrel_result <- readRDS(paste(file_prefix, "PCRelate", "RDS", sep="."))
grm <- pcrelateToMatrix(pcrel_result)

chunkSize <- length(all_covar_combinations) %/% clusterSize(cl)

cat(nrow(scan_annot), "rows in scan annotation df\n")
cat(length(sample.id), "sample names\n")
cat(nrow(grm), "rows in loaded GRM\n")

ptm <- proc.time()
models <- foreach (covars=all_covar_combinations,
         .packages=c("GENESIS", "stringr"),
         .errorhandling="stop",
         .options.mpi=list(chunkSize=chunkSize)) %dopar% {
        tryCatch({
            model <- fitNullModel(scan_annot, outcome=endpoint,
                                  covars=covars,
                                  sample.id=sample.id, cov.mat=grm, family=family)
            model$covars <- covars
            model
        }, error=function(cnd) list(error=cnd$message, converged=FALSE))
    }
exec_time <- proc.time() - ptm
models %>% compact("error") %>% map_chr("error") -> error_messages
models %>% keep(~ is.null(.x$error)) -> models
models %>% discard("converged") -> failed_models
models %>% keep("converged") -> succeeded_models
cat(length(all_covar_combinations), "covar combinations tested in", exec_time[["elapsed"]], "seconds\n")
cat(length(error_messages), "died with an error;", length(failed_models), "failed to converge;", length(succeeded_models), "succeeded\n")
cat("error summary:\n")
print(enframe(error_messages) %>% count(value))
if (length(succeeded_models) == 0) {
  stop("no models succeeded")
}
map_dbl(succeeded_models, "AIC") %>% compact %>% which.min %>% pluck(succeeded_models, .) -> model
cat("Best combination among succeeded models:", model$covars, "\n")

saveRDS(model, paste(file_prefix, paste0(endpoint, subset_tag), "null", "RDS", sep="."))

closeCluster(cl)
mpi.quit()
NULL

