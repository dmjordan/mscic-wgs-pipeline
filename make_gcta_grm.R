library(tibble)
library(forcats)
library(magrittr)
library(dplyr)
library(readr)

pcrelateToTibble <- function (infile, scaleKin=2) {
    pcrel <- readRDS(infile)
    pcrel$kinSelf %>% transmute(ID1=as_factor(ID),
                                ID2=as_factor(ID),
                                nsnp=nsnp,
                                kin=0.5*(1+f)) -> kinSelf
    pcrel$kinBtwn %>% transmute(ID1=factor(ID1, levels=levels(kinSelf$ID1)),
                                ID2=factor(ID2, levels=levels(kinSelf$ID2)),
                                nsnp=nsnp,
                                kin=kin) -> kinBtwn
    bind_rows(kinSelf, kinBtwn) %>% mutate(kin=scaleKin*kin)
}

writeGctaGrm <- function(grm_tbl, out_prefix) {
    grm_tbl %>%
      mutate(ID1=as.integer(ID1),
             ID2=as.integer(ID2),
             kin=formatC(kin, digits=4, format="f")) %>%
      arrange(ID1, ID2) %>%
      write_tsv(paste(out_prefix, "grm", "gz", sep="."), col_names=FALSE)

    tibble(FID=levels(grm_tbl$ID1), IID=levels(grm_tbl$ID1)) %>%
      write_tsv(paste(out_prefix, "grm", "id", sep="."), col_names=FALSE)
}

pcrelateToTibble(snakemake@input[[1]]) %>% writeGctaGrm(snakemake@wildcards[["prefix"]])
