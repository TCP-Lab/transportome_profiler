library(tidyverse)

## ------ Functions ----------> 
requireNamespace("DESeq2")

ccat <- function(...) {
  cat(paste0(...))
}


clean_data <- function(data, col = "Ensembl_ID") {
  
  give_first_split <- function(x, string) {
    res <- str_split_1(x, string)[1]
    
    res
  }
  
  wrapper <- function(x) {
    n <- sapply(x[[col]], give_first_split, string = "\\.")
    names(n) <- NULL
    x[[col]] <- n
    
    x
  }
  
  data <- lapply(data, wrapper)
  
  names(data) <- sapply(names(data), give_first_split, "\\.")
  
  data
}

split_gtex_data <- function(gtex_data, gtex_metadata) {
  gtypes <- unique(unlist(gtex_metadata[,"_primary_site"]))
  
  res <- list()
  for (gtype in gtypes) {
    samples <- gtex_metadata$Sample[gtex_metadata[,"_primary_site"] == gtype]
    
    ll <- length(samples)
    samples <- samples[samples %in% colnames(gtex_data)]
    ccat("Found ", length(samples), " for type ", gtype, ". Was ", ll, "\n")
    
    res[[gtype]] <- gtex_data[, c("sample", samples)]
    
  }
  
  tot_len <- sum(unlist(map(res, ncol))) - length(res)
  
  ccat(
    "Split ", tot_len, " entries of ", ncol(gtex_data) - 1,
    " total (", round(tot_len / (ncol(gtex_data) - 1) *100, 3), "%)")
  
  res
}

merge_on_rownames <- function(left, right, all.left = FALSE, all.right = FALSE, verbose = TRUE) {
  ll <- nrow(left)
  lr <- nrow(right)
  
  res <- merge(left, right, all.x = all.left, all.y = all.right, by = "row.names", sort = FALSE)
  row.names(res) <- res$Row.names
  res$Row.names <- NULL
  
  ccat("Merged left ", ll, " right ", lr, " >> ", nrow(res), "\n")
  
  res
}

run_deseq_on_frames <- function(healthy, tumor) {
  # Make a metadata dataframe
  meta <- data.frame(
    sample = c(colnames(healthy), colnames(tumor)),
    status = factor(c(rep("healthy", ncol(healthy)), rep("tumor", ncol(tumor))), levels = c("healthy", "tumor"))
  )
  
  meta |> column_to_rownames("sample") -> meta

  # Merge the two frames
  data <- merge_on_rownames(healthy, tumor)
  
  # Drop NAs
  data <- na.omit(data)
  
  # Make all values integers
  data <- round(data)
  
  # Sort the columns and rows to be right for DESeq2
  sorted_samples <- sort(colnames(data))
  data <- data[,sorted_samples]
  meta <- meta[sorted_samples,,drop=FALSE]

  assert(
    all(rownames(meta) == colnames(data)),
    "Different rownames and colnames between data and metadata"
  )
  
  data_matrix <- as.matrix(data)
  rownames(data_matrix) <- row.names(data)
  
  # Run DESeq2
  res <- DESeq2::DESeqDataSetFromMatrix(
    data_matrix, meta, ~ status
  ) |> DESeq2::DESeq()
  
  res
}

extract_result <- function(deseq_data, contrast) {
  df <- DESeq2::results(deseq_data, contrast = contrast, format = "DataFrame")
  
  df |> as.data.frame() -> df
  
  df
}

save_to_disk <- function(data, out_dir, name, rowname_col = "id") {
  name <- paste0(name, ".csv")
  ccat("Writing a ", nrow(data), " by ", ncol(data), " dataframe to ", file.path(out_dir, name), "\n" )
  
  data |> rownames_to_column(rowname_col) -> data
  write_csv(data, file.path(out_dir, name))
}

run_deseq <- function(tcga_data, gtex_data, tcga_gtex_matches, out_dir, skip_completed = TRUE) {
  present_files <- list.files(out_dir)
  for (i in seq_along(tcga_gtex_matches$TCGA)) {
    tcga_id <- tcga_gtex_matches$TCGA[i]
    gtex_id <- tcga_gtex_matches$GTEX[i]
    
    if (any(startsWith(present_files, tcga_id))) {
      ccat("Skipped generation of ", tcga_id, " vs gtex ", gtex_id, " as result was detected.\n")
      next
    }
    ccat("Runnig DESeQ2 on tcga: ", tcga_id, " and gtex: ", gtex_id, "\n")
    
    healthy <- gtex_data[[gtex_id]]
    tumor <- tcga_data[[tcga_id]]
    # Send the IDs to rownames
    healthy |> column_to_rownames("sample") -> healthy
    tumor |> column_to_rownames("Ensembl_ID") -> tumor
    
    deseq_obj <- run_deseq_on_frames(healthy, tumor)
    
    res <- extract_result(deseq_obj, c("status", "tumor", "healthy"))
    
    save_to_disk(res, out_dir, paste0(tcga_id, "_vs_", gtex_id))
  }
}

# <------ END Functions --------

# If you are running this from RStudio, you can skip this >>>>>>>>>>>>>>>>>>
if (sys.nframe() == 0L) {
  # Parsing arguments
  requireNamespace("argparser")

  parser <- argparser::arg_parser("Run DESeq2 on the data")

  parser |>
    argparser::add_argument(
      "expression_matrix", help="Expression matrix to run DESeq on.", type="character"
    ) |>
    argparser::add_argument(
      "output_directory", help="Folder to save the output directories in",
      type = "character"
    ) |>
    argparser::add_argument(
      "gtex_phenotype", help = "File with the GTEx phenotypes used to split the data with",
      type = "character"
    )|>
    argparser::add_argument(
      "tcga_metadata", help = "File with the TCGA metadata used to split the data with",
      type = "character"
    ) -> parser

  args <- argparser::parse_args(parser)
}
# To here <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<


MATCHES <- as.data.frame(rbind(
  c("LAML", "Blood"),
  c("ACC", "Adrenal Gland"),
  c("BLCA", "Bladder"),
  c("LGG", "Brain"),
  c("BRCA", "Breast"),
  c("CESC", "Cervix Uteri"),
  c("CHOL", "Liver"), # ??
  c("COAD", "Colon"),
  c("GBM", "Brain"),
  c("HNSC", "Skin"),
  c("KICH", "Kidney"),
  c("KIRC", "Kidney"),
  c("KIRP", "Kidney"),
  c("LIHC", "Liver"),
  c("LUAD", "Lung"),
  c("LUSC", "Lung"),
  c("DLBC", "Blood"),
  #c("MESO", ""), # ??
  c("OV", "Ovary"),
  c("PAAD", "Pancreas"),
  c("PCPG", "Kidney"),
  c("PRAD", "Prostate"),
  c("READ", "Colon"),
  #c("SARC", ""), # ??
  c("STAD", "Stomach"),
  c("TGCT", "Testis"),
  c("THYM", "Thyroid"),
  c("UCS", "Uterus"),
  c("UCEC", "Uterus")
  #c("UVM", "") # ??
))
colnames(MATCHES) <- c("TCGA", "GTEX")

all_data <- load_data("~/Files/data/tcga/TCGA/") 
data <- clean_data(all_data)
rm(all_data)

gtex_meta <- read_tsv("./GTEX_phenotype")
expression_matrix <- read_tsv(args$expression_matrix)
gtex_data_split <- split_gtex_data(expression_matrix, gtex_meta)
tcga_data_plit <- 
rm(gtex_data)
gtex_data_split <- clean_data(gtex_data_split, "sample")

run_deseq(data, gtex_data_split, MATCHES, "~/Files/data/tcga/TCGA_dea_output/")

