#!/bin/bash
## -*- mode:R -*-

## Walltime (minutes). Should be the time required without the benefit
## of caching.
#BSUB -W 720
## Parallelization happens via job submission, not multiprocessing
#BSUB -n 1
## Memory. The dream code has some unknown memory spikes.
#BSUB -M 50GB
#BSUB -R "rusage[mem=50GB]"

'\' >/dev/null 2>&1 || true
## This is bash code to set up the environment

module load R/4.0.3
module load python/3.7.3

## Bash setup code ends here
Rscript - "$@" <<"EOF";
invisible('#')
## R code starts after this line

suppressPackageStartupMessages({
    library(qs)
    library(treemap)
    library(R.matlab)
    library(goseq)
    library(topGO)
    library(org.Hs.eg.db)
    library(GO.db)
    library(Rgraphviz)
    library(rlang)
    library(readr)
    library(openxlsx)
    library(parallel)
    library(rms)
    library(broom)
    library(glmmLasso)
    library(tidyverse)
    filter <- dplyr::filter
    select <- dplyr::select
    rename <- dplyr::rename
    matches <- dplyr::matches
    library(magrittr)
    library(purrr)
    library(stringr)
    library(fs)
    path <- fs::path
    library(rex)
    library(assertthat)
    library(GenomicRanges)
    library(GenomicFeatures)
    genes <- GenomicFeatures::genes
    library(ensembldb)
    library(limma)
    library(edgeR)
    library(lme4)
    library(variancePartition)
    ## Make sure we get the S4 generic
    eBayes <- variancePartition::eBayes
    library(rctutils)
    library(SummarizedExperiment)
    library(splines)
    library(ggplot2)
    library(scales)
    library(treemapify)
    library(forcats)
    library(withr)
    library(BiocParallel)
    library(memoise)
    library(cachem)
    library(AnnotationHub)
    library(EnrichmentBrowser)
    library(topGO)
    library(ggrepel)
    library(retry)
    library(iterators)
    library(goseq)
    library(pryr)
    library(processx)
    run <- processx::run
    library(future)
    library(batchtools)
    library(future.batchtools)
    library(future.apply)
    library(BiocParallel.FutureParam)
    plan(multicore, workers = 10)
    register(MulticoreParam(30))
})

## Make sure we're on R 4.x
assert_that(compareVersion(R.version$major, "4") >= 0)

## Make sure memoise has the feature we need
assert_that("omit_args" %in% names(formals(args(memoise::memoise))))

options(
    # Export objects up to 2 GB
    future.globals.maxSize= 3*1024^3,
    # Let future fork inside RStudio
    future.fork.enable = TRUE
)

workdir <- path("/sc/arion/projects/mscic1/results/Ryan/RNAseq_analysis")
setwd(workdir)
plotdir <- path(workdir, "plots")
dir_create(plotdir)
datadir <- path(workdir, "saved_data")
dir_create(datadir)
results_dir <- path(workdir, "results")
dir_create(results_dir)

## source(path(workdir, "rnaseq-utils.R"))
## source(path(workdir, "force_assign.R"))

## Fill in NAs in first vector with non-NA values from subsequent
## vectors
na_fill <- function(x, ...) {
    if (length(x) == 0) {
        return(x)
    }
    fillers <- list(...)
    for (filler in fillers) {
        if (! length(filler) %in% c(1, length(x))) {
            warning("Length of filler is different from length of x")
        }
        filler <- rep(filler, length.out = length(x))
        to_fill <- is.na(x)
        if (any(to_fill)) {
            x[to_fill] <- filler[to_fill]
        } else {
            return(x)
        }
    }
    return(x)
}

## For use inside rex, because capture_group(n) doesn't work with
## stringr
capture_backref <- function(n) {
    rex:::p("\\", n)
}

getgo <- memoise(goseq::getgo)

## ## Force topGO::getSigGroups to use BiocParallel
## replace_in_package("topGO", ".sigGroups.classic", function (GOlist, test.stat)
## {
##     sigList <- bplapply(GOlist, function(term) {
##         members(test.stat) <- term
##         return(runTest(test.stat))
##     }, BPPARAM = BiocParallel::bpparam())
##     assertthat::assert_that(length(sigList) > 0)
##     sigList <- simplify2array(sigList, higher = (simplify == "array"))
##     return(sigList)
## })

## ## Memoise some internal topGO functions
## for (i in c("buildGOgraph.topology", "buildLevels", "mapGenes2GOgraph")) local({
##     f <- get(i, envir = ns_env("topGO"))
##     if (!is.memoised(f)) {
##         memo_f <- memoise(f)
##         replace_in_package("topGO", i, memo_f)
##     }
## })

## ## Register BPPARAM, evaluate expr, then undo registering BPPARAM
## with_bpparam <- function(BPPARAM, expr) {
##     bpreg <- environment(BiocParallel:::.register)$.self
##     current_registered <- bpreg$bpparams
##     withr::defer({
##         bpreg$bpparams <- current_registered
##     })
##     register(BPPARAM)
##     return(expr)
## }

## geneSig: vector of TRUE/FALSE significance, named by Ensembl gene
## ID. Names are the gene universe.
run_topGO_Fisher <- function(geneSig,
                             ontologies = c("BP", "MF", "CC"),
                             adjust_method="BH",
                             getgo_genome = "hg", getgo_id = "ensGene",
                             ## Currently does nothing
                             BPPARAM = BiocParallel::bpparam())
{
    ## Convert to 1/0 for Sig/NS
    geneSig <- setNames(as.numeric(as.logical(geneSig)), names(geneSig))
    assert_that(
        length(geneSig) > 1,
        !any(is.na(geneSig)),
        is_named(geneSig)
    )
    ontologies <- toupper(ontologies)
    assert_that(all(ontologies %in% c("BP", "MF", "CC")))
    ens_gene_map <- getgo(
        genes = names(geneSig),
        genome = getgo_genome,
        id = getgo_id,
        fetch.cats = str_c("GO:", toupper(ontologies))
    )
    test_stat = new(
        "classicCount", testStatistic = GOFisherTest,
        name = "Fisher test"
    )
    tgdata_list <- lapply(ontologies, function(onto) {
        new("topGOdata", description = paste0("Enrichment in ", onto),
            ontology = onto,
            allGenes = geneSig,
            geneSel = function(x) x != 0,
            nodeSize = 10,
            annotationFun = annFUN.gene2GO,
            gene2GO = ens_gene_map)
    }## , BPPARAM = BPPARAM
    )
    names(tgdata_list) <- ontologies
    node_counts <- vapply(tgdata_list, function(x) length(nodes(graph(x))), numeric(1))
    ## with_bpparam(BPPARAM, {
        res_list <- lapply(tgdata_list, getSigGroups, test.stat = test_stat)
    ## })
    res_tables <- mapply(GenTable, object = tgdata_list, classic = res_list, topNodes = node_counts, SIMPLIFY = FALSE)
    res_full_table <- res_tables %>%
        bind_rows(.id = "Ontology") %>%
        as_tibble %>%
        mutate(
            ## Some p-values are "< 1e-30" or similar, which we
            ## replace with just the number, as an upper bound.
            P.Value = classic %>%
                str_replace(rex(start, "<", maybe(spaces)), "") %>%
                as.numeric,
            classic = NULL,
            adj.P.Val = p.adjust(P.Value, method = adjust_method),
            Fold_Enrichment = Significant / Expected,
            N = n()
        )
    res_full_table
}
## run_topGO_Fisher <- memoise(run_topGO_Fisher_internal)

run_gofigure <- function(res_table,
                         ## Optional: default is to run in a temp dir
                         ## that gets deleted after (also deletes the
                         ## plots)
                         outdir,
                         outfile_suffix = str_c("GOFigure_", similarity_cutoff),
                         similarity_cutoff = 0.5,
                         adj.P.Threshold = 0.05,
                         gofigure_path = "/sc/arion/projects/mscic1/data/GOFigure/GO-Figure/gofigure.py"
                         ) {
    assert_that(file_exists(gofigure_path))
    if (missing(outdir)) {
        outdir <- tempfile()
        withr::defer({
            if (dir_exists(outdir)) {
                dir_delete(outdir)
            }
        })
        dir.create(outdir)
    }
    assert_that(dir_exists(outdir))
    sig_res_table <- res_table %>%
        mutate(adj.P.Val = pmax(adj.P.Val, 10^-300, na.rm = TRUE)) %>%
        filter(adj.P.Val <= adj.P.Threshold)
    ## Put in "standard input" format:
    ## https://gitlab.com/evogenlab/GO-Figure#standard-input
    standard_input_table <- sig_res_table %>%
        select(GO.ID, adj.P.Val)
    resList <- list()
    for (onto in c("biological_process", "cellular_component", "molecular_function")) {
        resList[[onto]] <- NULL
    }
    if (nrow(sig_res_table) == 0) {
        warning("No significant GO terms found")
        return(resList)
    }
    with_tempfile("tempf", {
        write_tsv(standard_input_table, tempf, col_names = FALSE)
        gofigure_cmd <- c(
            "python3",
            gofigure_path,
            "--input", tempf,
            "--input_type", "standard",
            "--output", outdir,
            "--ontology", "all",
            "--file_type", "pdf",
            "--outfile_appendix", outfile_suffix,
            "--similarity_cutoff", similarity_cutoff
        )
        result <- processx::run(
            command = gofigure_cmd[1],
            args = gofigure_cmd[-1]
        )
    })
    resList <- list()
    for (onto in c("biological_process", "cellular_component", "molecular_function")) {
        fname <- path(outdir, str_c(onto, "_full_table_", outfile_suffix, ".tsv"))
        if (file_exists(fname)) {
            resList[[onto]] <- suppressMessages(read_tsv(fname))
        } else {
            warning("Not enough significant GO terms for ", onto)
        }
    }
    resList
}

## Function for making an empty "placeholder" plot with no content.
empty_plot <- function(...) {
    plot(0,type='n',axes=TRUE, xaxt='n', yaxt = "n", xlab = "", ylab = "", ...)
}

## table1 and table2 are GO-Figure results tables with lower and
## higher similarity indices, respectively.
plot_gofigure_module_treemap <- function(table1, table2, title = "GO-Figure tree map", ...) {
    if (!is.data.frame(table1) && is.list(table1)) {
        table1 <- bind_rows(table1)
    }
    assert_that(is.data.frame(table1))
    if (!is.data.frame(table2) && is.list(table2)) {
        table2 <- bind_rows(table2)
    }
    assert_that(is.data.frame(table2))
    table2 %<>% mutate(
        eliminated = ! `Cluster member` %in% `Cluster representative`,
        real_representative = table1$`Cluster representative`[match(`Cluster member`, table1$`Cluster member`)]
    )
    clusters <- table2 %>%
        filter(!eliminated) %>%
        mutate(
            cluster_name = `Cluster member description`[match(real_representative, `Cluster member`)] %>%
                na_fill(""),
                ## factor(levels = `Cluster member`) %>% as.character %>% na_fill(""),
            log10_p.value = pmax(0, -log10(as.numeric(`Member P-value`)), na.rm = TRUE)
        )
    ## Note: need to implement subgroups for ggplot code
    ## p <- ggplot(clusters) +
    ##     aes(
    ##         area = log10_p.value,
    ##         fill = cluster_name,
    ##         label = cluster_name
    ##     ) +
    ##     geom_treemap(
    ##         start = "topleft"
    ##     ) +
    ##     guides(fill = "none") +
    ##     geom_treemap_text(
    ##         start = "topleft",
    ##         colour = "black",
    ##         place = "centre",
    ##         reflow = TRUE
    ##     ) +
    ##     ggtitle(title)
    if (nrow(table1) > 0 && nrow(table2) > 0) {
        treemap(
            clusters,
            index = c("cluster_name","Cluster member description"),
            vSize = "log10_p.value",
            type = "categorical",
            vColor = "cluster_name",
            title = title,
            inflate.labels = FALSE,      # set this to TRUE for space-filling group labels - good for posters
            lowerbound.cex.labels = 0,   # try to draw as many labels as possible (still, some small squares may not get a label)
            bg.labels = "#CCCCCCAA",     # define background color of group labels
            # "#CCCCCC00" is fully transparent, "#CCCCCCAA" is semi-transparent grey, NA is opaque
            position.legend = "none",
            ...
        )
    } else {
        stop("No significant clusters. Cannot plot tree map")
    }
}

## system.time(full_de_result_table <- data.table::fread(path(datadir, "full_de_result_table.csv")))
{
    message("Loading DE gene lists")
    de_result_table_list <- qread(path(datadir, "de_result_table_list.qs"), nthreads = data.table::getDTthreads())
    all_contrasts_table <- de_result_table_list %>%
        lapply(. %>% .[1,c("Source", "Name", "ModelName", "CellAdjust", "CellInteraction", "Contrast")]) %>%
        bind_rows(.id = "TableName")
    plan(multicore, workers = data.table::getDTthreads())
    geneSig_direction_list <- future_lapply(
        de_result_table_list,
        . %>%  group_by(gene_id) %>%
        summarise(Significant_Direction = any(adj.P.Val <= 0.05) * sign(logFC)) %>%
        deframe
    )
    rm(de_result_table_list)
    invisible(gc())
}

{
    message("Generating directional and non-directional gene lists")
    geneSig_lists <- list(
        Nondirectional = lapply(geneSig_direction_list, function(x) x != 0),
        Up = lapply(geneSig_direction_list, function(x) x > 0),
        Down = lapply(geneSig_direction_list, function(x) x < 0)
    )
}

topGO_res_lists <- list()
for (i in names(geneSig_lists)) {
    message("Running topGO for direction ", i)
    ## Construct an empty list with the right length and names
    topGO_res_lists[[i]] <- rep(list(NULL), length = length(geneSig_lists[[i]]))
    names(topGO_res_lists[[i]]) <- names(geneSig_lists[[i]])
    message("Checking for cached topGO results")
    topGO_datadir <- path(datadir, "topGO")
    topGO_data_files <- path(
        topGO_datadir,
        str_c("topGO-", names(geneSig_lists[[i]]), "-" ,i, ".RDS")
    ) %>%
        set_names(names(geneSig_lists[[i]]))
    dir_create(unique(path_dir(topGO_data_files)))
    already_ran <- file_exists(topGO_data_files)
    if (any(already_ran)) {
        message("Loading ", sum(already_ran), " previously run topGO results.")
        plan(multicore, workers = data.table::getDTthreads())
        topGO_res_lists[[i]][already_ran] <- future_lapply(
            topGO_data_files[already_ran],
            readRDS,
            future.chunk.size = 20
        )
    }
    need_to_run <- sapply(topGO_res_lists[[i]], function(x) is.null(x) || is(x, "try-error")) %>%
        which %>% names
    if (length(need_to_run) > 0) {
        message("Need to run topGO on ", length(need_to_run), " tables")
        ## Now actually run the ones that haven't been memoised
        ## plan(multicore, workers = 30)
        chunk_size = min(10, length(need_to_run)/10)
        plan_cluster_jobs <- tweak(
            batchtools_custom,
            resources = list(
                walltime = 60 * (chunk_size * 3 + 5),
                threads_per_job = 1,
                memory = "3GB",
                queue = "express"
            ),
            workers = ceiling(length(need_to_run) / chunk_size)
        )
        plan(plan_cluster_jobs)
        topGO_res_lists[[i]][need_to_run] <- future_mapply(
            function(gs, fname) {
                res <- run_topGO_Fisher(gs)
                dir_create(path_dir(fname))
                saveRDS(res, fname)
                res
            },
            gs = geneSig_lists[[i]][need_to_run],
            fname = topGO_data_files[need_to_run],
            future.chunk.size = chunk_size,
            SIMPLIFY = FALSE
        )
        need_to_run <- sapply(topGO_res_lists[[i]], function(x) is.null(x) || is(x, "try-error")) %>%
            which %>% names
        assert_that(!any(need_to_run))
    }
}

gofigure_0.1_res_lists <- list()
gofigure_0.7_res_lists <- list()
for (i in names(topGO_res_lists)) {
    message("Running GO-Figure for direction ", i, " (or loading previously run results)")
    chunk_size = min(100, length(topGO_res_lists[[i]])/30)
    ## plan_cluster_jobs <- tweak(
    ##     batchtools_custom,
    ##     resources = list(
    ##         walltime = 60 * (chunk_size * 0.3 + 5),
    ##         threads_per_job = 1,
    ##         memory = "3GB",
    ##         queue = "express"
    ##     ),
    ##     workers = ceiling(length(topGO_res_lists) / chunk_size)
    ## )
    ## plan(plan_cluster_jobs)
    plan(multicore, workers = data.table::getDTthreads())
    gofigure_0.1_data_file <- path(datadir, "gofigure", str_c("gofigure_0.1_", i, "_res_list.RDS"))
    gofigure_0.1_res_lists[[i]] <- tryCatch({
        x <- readRDS(gofigure_0.1_data_file)
        assert_that(length(x) == length(topGO_res_lists[[i]]))
        message("Loaded previously-computed GOfigure results from file: ", gofigure_0.1_data_file)
        x
    }, error = function(e) {
        message("Running GOfigure with si=0.1 for ", i)
        res_list <- future_lapply(
            topGO_res_lists[[i]],
            function(...) try(run_gofigure(...)),
            similarity_cutoff = 0.1,
            future.seed = TRUE,
            future.chunk.size = 30 # chunk_size
        )
        ## Only save the data file if all succeeded
        succeeded <- !sapply(res_list, function(x) is.null(x) || is(x, "try-error"))
        if (all(succeeded)) {
            message("Saving GOfigure data file: ", gofigure_0.1_data_file)
            dir_create(path_dir(gofigure_0.1_data_file))
            saveRDS(res_list, gofigure_0.1_data_file)
        } else {
            message("Not saving GOfigure data file because some runs had errors: ", gofigure_0.1_data_file)
        }
        res_list
    })
    gofigure_0.7_data_file <- path(datadir, "gofigure", str_c("gofigure_0.7_", i, "_res_list.RDS"))
    gofigure_0.7_res_lists[[i]] <- tryCatch({
        x <- readRDS(gofigure_0.7_data_file)
        assert_that(length(x) == length(topGO_res_lists[[i]]))
        message("Loaded previously-computed GOfigure results from file: ", gofigure_0.7_data_file)
        x
    }, error = function(e) {
        message("Running GOfigure with si=0.7 for ", i)
        res_list <- future_lapply(
            topGO_res_lists[[i]],
            function(...) try(run_gofigure(...)),
            similarity_cutoff = 0.7,
            future.seed = TRUE,
            future.chunk.size = 30 # chunk_size
        )
        succeeded <- !sapply(res_list, function(x) is.null(x) || is(x, "try-error"))
        if (all(succeeded)) {
            message("Saving GOfigure data file: ", gofigure_0.7_data_file)
            dir_create(path_dir(gofigure_0.7_data_file))
            saveRDS(res_list, gofigure_0.7_data_file)
        } else {
            message("Not saving GOfigure data file because some runs had errors: ", gofigure_0.7_data_file)
        }
        res_list
    })
}

## For paths, do one folder per model, one PDF per contrast with pages
## for each cell type variant. For transitions, do a 2-level
## hierarchy, with transition model and then day number.
for (direction in names(gofigure_0.1_res_lists)) {
    message("Creating GO-Figure tree map plots for direction ", direction)
    gofigure_plot_file_table <- all_contrasts_table %>%
        mutate(
            OutDir = path(
                plotdir, "GOfigure",
                ModelName %>%
                str_replace(
                    rex(capture("transitionName", anything), "_day"),
                    "\\1/day"
                ) %>%
                str_replace(
                    rex("_" %if_prev_is% "Post_COVID19"),
                    "/"
                )
            ),
            OutFile = path(OutDir, str_c(Contrast, "_", direction, ".pdf")),
            Plot_Title = str_c(
                ModelName,
                case_when(
                    !is.na(CellInteraction) ~ str_c(".", CellInteraction),
                    CellAdjust ~ ".CellAdjust",
                    TRUE ~ ""
                ),
                if_else(
                    ModelName == Contrast,
                    "",
                    Contrast
                ),
                " ", direction
            ),
            ## Plot_Title = str_c(TableName, " ", direction),
            table1 = gofigure_0.1_res_lists[[direction]][TableName] %>% lapply(bind_rows),
            table2 = gofigure_0.7_res_lists[[direction]][TableName] %>% lapply(bind_rows)
        ) %>%
        arrange(OutFile, CellAdjust, !is.na(CellInteraction), CellInteraction)
    dir_create(unique(gofigure_plot_file_table$OutDir))
    gofigure_plot_file_table_list <- split(gofigure_plot_file_table, gofigure_plot_file_table$OutFile)
    plan(multicore, workers = data.table::getDTthreads())
    future_lapply(gofigure_plot_file_table_list, function(x) {
        fname <- unique(x$OutFile)
        assert_that(is_string(fname))
        message("Creating GOfigure plots in ", deparse1(fname))
        with_pdf(fname, {
            for (i in seq_len(nrow(x))) {
                ## Do the plot
                table_name <- x$TableName[i]
                plot_title <- x$Plot_Title[i]
                table1 <- x$table1[[i]]
                table2 <- x$table2[[i]]
                if (nrow(table1) > 0 && nrow(table2) > 0) {
                    plot_gofigure_module_treemap(table1, table2, plot_title)
                } else {
                    ## If nothing is significant, make an empty plot,
                    ## since we still need to generate a PDF page.
                    empty_plot(
                        main = plot_title,
                        sub = "[No significant clusters]",
                        font.main = 1,
                        font.sub = 1
                    )
                }
            }
        })
    })
}

## R code ends here
EOF <- NULL
EOF
