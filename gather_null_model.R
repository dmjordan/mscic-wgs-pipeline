library(tidyverse)

models <- lapply(snakemake@input, readRDS)

models %>% compact("error") %>% map_chr("error") -> error_messages
models %>% keep(~ is.null(.x$error)) -> models
models %>% discard(~ifelse(is.null(.x$converged), FALSE, .x$converged)) -> failed_models
models %>% keep(~ifelse(is.null(.x$converged), FALSE, .x$converged)) -> succeeded_models
cat("of", length(error_messages) + length(models), "models:", length(error_messages), "died with an error;", length(failed_models), "failed to converge;", length(succeeded_models), "succeeded\n")
cat("error summary:\n")
print(enframe(error_messages) %>% count(value))
map_dbl(succeeded_models, "AIC") %>% compact %>% which.min %>% pluck(succeeded_models, .) -> model
cat("Best combination among succeeded models:", model$covars, "\n")

saveRDS(model, snakemake@output[[1]])
