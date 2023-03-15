library(tidyverse)

requireNamespace("DESeq2")

load_data <- function(input_folder) {
  files <- list.files(input_folder)
  file_paths <- file.path(input_folder, files)
  
  res <- list()
  for (i in seq_along(files)) {
    res[[files[i]]] <- read_tsv(file_paths[i], show_col_types = FALSE)
  }
  
  res
}

all_data <- load_data("~/Files/data/tcga/TCGA/") 

MATCHES <- rbind(
  c("TCGA-LAML", "Blood"),
  c("TCGA-ACC", "Adrenal Gland"),
  c("TCGA-BLCA", "Bladder"),
  c("TCGA-LGG", "Brain"),
  c("TCGA-BRCA", "Breast"),
  c("TCGA-CESC", "Cervix Uteri"),
  c("TCGA-CHOL", "Liver"), # ??
  c("TCGA-LCML", "Blood"),
  c("TCGA-COAD", "Colon"),
  c("TCGA-ESCA", "Esophagus"),
  c("TCGA-GBM", "Brain"),
  c("TCGA-HNSC", "Skin"),
  c("TCGA-KICH", "Kidney"),
  c("TCGA-KIRC", "Kidney"),
  c("TCGA-KIRP", "Kidney"),
  c("TCGA-LIHC", "Liver"),
  c("TCGA-LUAD", "Lung"),
  c("TCGA-LUSC", "Lung"),
  c("TCGA-DLBC", "Blood"),
  #c("TCGA-MESO", ""), # ??
  c("TCGA-OV", "Ovary"),
  c("TCGA-PAAD", "Pancreas"),
  c("TCGA-PCPG", "Kidney"),
  c("TCGA-PRAD", "Prostate"),
  c("TCGA-READ", "Colon"),
  #c("TCGA-SARC", ""), # ??
  c("TCGA-STAD", "Stomach"),
  c("TCGA-TGCT", "Testis"),
  c("TCGA-THYM", "Thyroid"),
  c("TCGA-UCS", "Uterus"),
  c("TCGA-UCEC", "Uterus"),
  #c("TCGA-UVM", "") # ??
)

clean_data <- function(data) {
  
  give_first_split <- function(x, string) {
    res <- str_split_1(x, string)[1]
    
    res
  }
  
  wrapper <- function(x) {
    n <- sapply(x$Ensembl_ID, give_first_split, string = "\\.")
    names(n) <- NULL
    x$Ensembl_ID <- n
    
    x
  }
  
  data <- lapply(data, wrapper)
  
  names(data) <- sapply(names(data), give_first_split, "\\.")
  
  data
}

data <- clean_data(all_data)

rm(all_data)

run_deseq <- function(data) {
  
}


# --------------

gtex_meta <- read_tsv("./GTEX_phenotype")

table(gtex_meta$`_primary_site`)

