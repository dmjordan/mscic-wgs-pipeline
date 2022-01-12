library(tidyverse)

models <- lapply(snakemake@input, readRDS)

models %>% compact("error") %>% map_chr("error") -> error_messages
models %>% keep(~ is.null(.x$error)) -> models
models %>% discard("converged") -> failed_models
models %>% keep("converged") -> succeeded_models
cat(length(all_covar_combinations), "covar combinations tested in", exec_time[["elapsed"]], "seconds\n")
cat(length(error_messages), "died with an error;", length(failed_models), "failed to converge;", length(succeeded_models), "succeeded\n")
cat("error summary:\n")
print(enframe(error_messages) %>% count(value))
map_dbl(succeeded_models, "AIC") %>% compact %>% which.min %>% pluck(succeeded_models, .) -> model
cat("Best combination among succeeded models:", model$covars, "\n")

saveRDS(model, snakemake@output[[1]])