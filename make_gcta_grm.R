library(tibble)
library(forcats)
library(magrittr)
library(dplyr)
library(readr)
library(GENESIS)

if (!require(genio)) { 
    install.packages("genio", repos="https://repo.miserver.it.umich.edu/cran/")
    library(genio)
}

pcrelateToNsnpsMatrix <- function (pcrel) {
    pcrel$kinSelf %>% transmute(ID1=ID,
                                ID2=ID,
                                value=nsnp) -> nSelf
    pcrel$kinBtwn %>% transmute(ID1=ID1,
                                ID2=ID2,
                                value=nsnp) -> nBtwn
    bind_rows(nSelf, nBtwn) %>% makeSparseMatrix
}

pcrelate <- readRDS(snakemake@input[[1]])
grm <- pcrelateToMatrix(pcrelate)
nsnps_matrix <- pcrelateToNsnpsMatrix(pcrelate)
ids <- tibble(fam=rownames(grm), id=rownames(grm))

write_grm(snakemake@wildcards[["prefix"]], 
          as.matrix(grm), 
          as.matrix(nsnps_matrix),
          ids)
          
