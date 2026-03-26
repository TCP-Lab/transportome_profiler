#' This script compares transcriptomics and proteomics using the
#' CPTAC datasets.

library(tidyverse)
library(progress)

TUMOR_TYPES <- c("BRCA", "CCRCC", "COAD", "GBM", "HNSCC", "LSCC", "LUAD", "OV", "PDAC", "UCEC")
ENSG_DATA <- read_csv(file.path("data", "in", "ensg_data.csv"))
CODING_GENES <- ENSG_DATA |> filter(gene_biotype == "protein_coding")
GENESETS <- jsonlite::read_json(file.path("data", "genesets.json"))

get_geneset_named <- function(data, x) {
    for (item in data) {
        if (item$name == x) {
            return(item)
        }
    }
    stop("Item not found")
}

TRANSPORTOME <- get_geneset_named(GENESETS, "whole_transportome")$data |> unlist()
CHANNELS <- get_geneset_named(GENESETS, "channels")$data |> unlist()
TRANSPORTERS <- get_geneset_named(GENESETS, "transporters")$data |> unlist()

load_data <- function(root_path) {
    data <- list()
    for (ttype in TUMOR_TYPES) {
        this_tt <- list()
        this_tt$tumor <- list()
        this_tt$normal <- list()
        # This is slightly duplicated but it's just too hard to make it concise
        cat(paste0("Loading ", ttype, " tumor proteomics data\n"))
        tryCatch({
            this_tt$tumor$proteomics <- read_tsv(file.path(root_path, paste0(ttype, "_tumor_proteomics.tsv")), show_col_types = FALSE)
        }, error = function(e) {cat(paste0("Failed to load ", ttype, " tumor proteomics data\n"))})
        
        cat(paste0("Loading ", ttype, " tumor rnaseq data\n"))
        tryCatch({
            this_tt$tumor$rnaseq <- read_tsv(file.path(root_path, paste0(ttype, "_tumor_RSEM.tsv")), show_col_types = FALSE)
        }, error = function(e) {cat(paste0("Failed to load ", ttype, " tumor rnaseq data\n"))})
        
        cat(paste0("Loading ", ttype, " normal proteomics data\n"))
        tryCatch({
            this_tt$normal$proteomics <- read_tsv(file.path(root_path, paste0(ttype, "_normal_proteomics.tsv")), show_col_types = FALSE)
        }, error = function(e) {cat(paste0("Failed to load ", ttype, " normal proteomics data\n"))})
        
        cat(paste0("Loading ", ttype, " normal rnaseq data\n"))
        tryCatch({
            this_tt$normal$rnaseq <- read_tsv(file.path(root_path, paste0(ttype, "_normal_RSEM.tsv")), show_col_types = FALSE)
        }, error = function(e) {cat(paste0("Failed to load ", ttype, " normal rnaseq data\n"))})
        
        data[[ttype]] <- this_tt
    }
    
    data
}

data <- load_data("./data/in/proteomics")

# Utility to detect is some data or other is missing or present
has <- function(data, type = c("tumor", "normal"), omic = c("proteomics", "transcriptomics")) {
    layer <- data[[type]]
    if (startsWith(omic, "p")) {
        return(!is.null(layer$proteomics))
    } else if (startsWith(omic, "t")) {
        return(!is.null(layer$rnaseq))
    } else {
        stop("Cannot find specified omic")
    }
}

has_tumor_seq <- partial(has, type = "tumor", omic = "trans")
has_tumor_prot <- partial(has, type = "tumor", omic = "prot")
has_normal_seq <- partial(has, type = "normal", omic = "trans")
has_normal_prot <- partial(has, type = "normal", omic = "prot")

# Utility to check if two frames have the same columns
check_samples <- function(this, that, verbose = FALSE) {
    out <- FALSE
    one_names <- names(this)
    two_names <- names(that)
    
    res <- c()
    
    one_ln <- length(one_names)
    two_ln <- length(two_names)
    res <- c(res, "First DF has", one_ln, "cols. The other has", two_ln, ".\n")
    
    inters <- intersect(one_names, two_names)
    one_uni <- setdiff(one_names, inters)
    two_uni <- setdiff(two_names, inters)
    
    if (length(inters) > 0 & length(one_uni) == 0 & length(two_uni) == 0) {
        res <- c(res, "All columns intersect. There are no unique columns.")
        out <- TRUE
    } else if (!verbose) {
        res <- c(res, "There are", length(inters), "intersecting cols and", length(one_uni), ",", length(two_uni), "unique columns in the first and second df, respectively.\n")
    } else {
        res <- c(
            res, "There are", length(inters), "intersecting cols:", inters, "\n",
            "There are", length(one_uni), "unique cols in the first df:", one_uni, "\n",
            "There are", length(two_uni), "unique cols in the second df:", two_uni, "\n"
        )
    }
    
    cat(paste0(res))
    
    invisible(out)
}

# Filters two dataframes and returns them as a list, keeping only intersecting columns
force_homogeneity <- function(this, that, force = NULL) {
    if (!is.null(force)) {
        stopifnot("Forced column(s) is not present in first DF" = all(force %in% names(this)))
        stopifnot("Forced column(s) is not present in second DF" = all(force %in% names(that)))
    }
    
    intersecting_cols <- intersect(names(this), names(that))
    
    this <- this |> select(all_of(intersecting_cols))
    that <- that |> select(all_of(intersecting_cols))
    
    return(list(this, that))
}

strip_version <- function(x) {
    str_split_i(x, "\\.", 1)
}

collapse_duplicate_genes <- function(data, id_col = "idx", method = "drop") {
    dups <- data[[id_col]] |> duplicated()
    if (any(dups)) {
        cat(paste0("Dropping ", sum(dups), " duplicates...\n"))
    }
    if (method == "drop") {
        data <- data[!dups, ]
    } else {
        stop("Unrecognized method")
    }
    
    data
}

purge_gene_versions <- function(data, id_col = "idx") {
    data[[id_col]] <- strip_version(data[[id_col]])
    data
}

select_coding <- function(data, id_col = "idx") {
    presence <- CODING_GENES$ensembl_gene_id %in% data[[id_col]]
    if (!all(presence)) {
        warning(paste0("There are ", sum(!presence), " coding genes not present in the input."))
    }
    
    data |> filter(.data[[id_col]] %in% CODING_GENES$ensembl_gene_id)
}

intersect_genes <- function(this, that, id_col = "idx") {
    this_genes <- this[[id_col]]
    that_genes <- that[[id_col]]
    
    common <- intersect(this_genes, that_genes)
    
    if (length(common) < 0.7 * length(this_genes)) {
        warning(paste0("Intersection of first DF left less than ", round(length(common) / length(this_genes) * 100, 2)), "% genes")
    }
    if (length(common) < 0.7 * length(that_genes)) {
        warning(paste0("Intersection of second DF left less than ", round(length(common) / length(that_genes) * 100, 2)), "% genes")
    }
    
    this <- this[this[[id_col]] %in% common, ]
    that <- that[that[[id_col]] %in% common, ]
    
    return(list(this, that))
}

intersect_genes(
    data$BRCA$tumor$proteomics |> purge_gene_versions() |> select_coding() |> collapse_duplicate_genes(),
    data$BRCA$tumor$rnaseq  |> purge_gene_versions() |> select_coding() |> collapse_duplicate_genes()
)

# Calculates the correlations between all genes in two data frames
calculate_correlations <- function(this, that, id_col = "idx") {
    if (any(is.na(this))) {
        warning("First dataframe has some NAs. Weird things might happen.")
    }
    if (any(is.na(that))) {
        warning("Second dataframe has some NAs. Weird things might happen.")
    }
    
    genes <- intersect(this[[id_col]], that[[id_col]])

    results <- list()
    pb <- progress_bar$new(total = length(genes))
    for (gene in genes) {
        this_gene <- this |> filter(idx == gene) |> select(!{{ id_col }}) |> unlist()
        that_gene <- that |> filter(idx == gene) |> select(!{{ id_col }}) |> unlist()
        
        result <- list()
        
        result$corr <- suppressWarnings(cor.test(this_gene, that_gene, method = "spearman")$estimate)
        result$mean <- mean(c(this_gene, that_gene))
        
        pb$tick()
        
        results[[gene]] <- result
    }
    final <- bind_rows(results)
    final$gene <- genes
    
    final
}

process_pair <- function(this, that, id_col = "idx", intersect_samples = TRUE) {
    check_samples(this, that)
    this <- this |> purge_gene_versions(id_col=id_col) |> select_coding(id_col=id_col) |> collapse_duplicate_genes(id_col=id_col)
    that <- that |> purge_gene_versions(id_col=id_col) |> select_coding(id_col=id_col) |> collapse_duplicate_genes(id_col=id_col)
    
    res <- force_homogeneity(this, that)
    this <- res[[1]]
    that <- res[[2]]
    
    res <- intersect_genes(drop_na(this), drop_na(that), id_col=id_col)
    this <- res[[1]]
    that <- res[[2]]
    
    return(calculate_correlations(this, that, id_col=id_col))
}

subset_with <- function(data, selected_genes, id_col = "idx") {
    data[strip_version(data[[id_col]]) %in% selected_genes, ]
}

mmerge <- function(this, that, id_col = "idx") {
    merge(this, that, by = id_col, all = TRUE)
}

process_batch <- function(data) {
    # Some calls here are duplicated to be slightly more explicit, and in case
    # we need to edit only some cases.
    
    results <- list()
    ## TUMOR DATA
    if (has_tumor_prot(data) & has_tumor_seq(data)) {
        results$tumor <- list()
        cat("Processing - TUMORS\n")
        # All genes
        cat("Processing Tumors - all\n")
        results$tumor$all <- process_pair(data$tumor$proteomics, data$tumor$rnaseq)
        # Only whole transportome
        cat("Processing Tumors - whole transportome\n")
        results$tumor$whole_transportome <- process_pair(
            subset_with(data$tumor$proteomics, TRANSPORTOME),
            subset_with(data$tumor$rnaseq, TRANSPORTOME)
        )
        # Only channels
        cat("Processing Tumors - channels\n")
        results$tumor$channels <- process_pair(
            subset_with(data$tumor$proteomics, CHANNELS),
            subset_with(data$tumor$rnaseq, CHANNELS)
        )
        # Only transporters
        cat("Processing Tumors - transporters\n")
        results$tumor$transporters <- process_pair(
            subset_with(data$tumor$proteomics, TRANSPORTERS),
            subset_with(data$tumor$rnaseq, TRANSPORTERS)
        )
    }
    
    if (has_normal_prot(data) & has_normal_seq(data)) {
        results$normal <- list()
        # All genes
        cat("Processing Tumors - all\n")
        results$normal$all <- process_pair(data$normal$proteomics, data$normal$rnaseq)
        # Only whole transportome
        cat("Processing Tumors - whole transportome\n")
        results$normal$whole_transportome <- process_pair(
            subset_with(data$normal$proteomics, TRANSPORTOME),
            subset_with(data$normal$rnaseq, TRANSPORTOME)
        )
        # Only channels
        cat("Processing Tumors - channels\n")
        results$normal$channels <- process_pair(
            subset_with(data$normal$proteomics, CHANNELS),
            subset_with(data$normal$rnaseq, CHANNELS)
        )
        # Only transporters
        cat("Processing Tumors - transporters\n")
        results$normal$transporters <- process_pair(
            subset_with(data$normal$proteomics, TRANSPORTERS),
            subset_with(data$normal$rnaseq, TRANSPORTERS)
        )
    }
    
    if (all(c(has_normal_prot(data), has_normal_seq(data), has_tumor_prot(data), has_tumor_seq(data)))) {
        results$combined <- list()
        all_proteomics <- mmerge(data$normal$proteomics, data$tumor$proteomics)
        all_transcriptomics <- mmerge(data$normal$rnaseq, data$tumor$rnaseq)
        
        # All genes
        cat("Processing Combined - all\n")
        results$combined$all <- process_pair(all_proteomics, all_transcriptomics)
        # Only whole transportome
        cat("Processing Combined - whole transportome\n")
        results$combined$whole_transportome <- process_pair(
            subset_with(all_proteomics, TRANSPORTOME),
            subset_with(all_transcriptomics, TRANSPORTOME)
        )
        # Only channels
        cat("Processing Combined - channels\n")
        results$combined$channels <- process_pair(
            subset_with(all_proteomics, CHANNELS),
            subset_with(all_transcriptomics, CHANNELS)
        )
        # Only transporters
        cat("Processing Combined - transporters\n")
        results$combined$transporters <- process_pair(
            subset_with(all_proteomics, TRANSPORTERS),
            subset_with(all_transcriptomics, TRANSPORTERS)
        )
    }
    
    results
}

calc_all_correlations <- function() {
    results <- list()
    for (ttype in TUMOR_TYPES) {
        cat(paste0("Processing batches for type ", ttype, ".\n"))
        tryCatch({
            results[[ttype]] <- process_batch(data[[ttype]])
        }, error = function(e) {print(paste0("Failed to process ", ttype, " Error: ", str(e)))})
    }
    
    results
}

correlation_results <- calc_all_correlations()

## - Plotting -

prepare_plot_data <- function(corrs) {
    noerr <- partial(try, silent = TRUE)
    flat_res <- list()
    i <- 1
    # Add the various identifiers
    for (ttype in names(corrs)) {
        for (status in c("tumor", "normal", "combined")) {
            for (test in c("all", "whole_transportome", "channels", "transporters")) {
                noerr({
                    corrs[[ttype]][[status]][[test]]$tumor_type <- ttype
                    corrs[[ttype]][[status]][[test]]$status <- status
                    corrs[[ttype]][[status]][[test]]$test <- test
                    
                    flat_res[[i]] <- corrs[[ttype]][[status]][[test]]
                    i <- i + 1
                })
            }
        }
    }
    
    # Collapse to a single frame
    bind_rows(flat_res)
}

plot_data <- prepare_plot_data(correlation_results)

x <- ggplot(plot_data, aes(x = status, fill = test, y = corr)) +
    geom_boxplot() +
    facet_wrap(facets = ~ tumor_type) +
    theme(legend.position = "bottom") +
    ylab("Spearmann's Correlation") +
    xlab("Cohort") +
    scale_fill_discrete(name = "Geneset")

pdf(
    file = file.path("data", "out", "transcriptomics_proteomics_plot.pdf"),
    width = 16, height = 9
)
print(x)
dev.off()
