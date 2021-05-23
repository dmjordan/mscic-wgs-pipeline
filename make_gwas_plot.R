suppressPackageStartupMessages({
  library(qqman)
  library(tidyverse)
})

args <- commandArgs(trailingOnly=TRUE)
endpoint <- args[[1]]

assoc <- read_delim(paste(endpoint, "GENESIS", "assoc", "txt", sep="."), delim=" ")

png(paste(endpoint, "GENESIS", "qq", "png", sep="."))
qq(assoc$Score.pval)
dev.off()

png(paste(endpoint, "GENESIS", "manhattan", "png", sep="."))
mutate(assoc, chr_numeric=as.numeric(as.factor(chr))) %>% filter(chr_numeric < 25) %>%
    manhattan(chr="chr_numeric", bp="pos", p="Score.pval", snp="variant.id",
                chrlabs=c(1:22,"X","Y"), col=c("blue4","orange3"))#,ylim=c(0,10))
dev.off()
