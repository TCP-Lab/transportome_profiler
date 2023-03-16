
# For ease of loading - i'm going to hate myself later
setwd("~/Files/repos/TCGA_dea/src/deseq")

source("./deseq_runner_fns.R")

MATCHES <- as.data.frame(rbind(
  c("TCGA-LAML", "Blood"),
  c("TCGA-ACC", "Adrenal Gland"),
  c("TCGA-BLCA", "Bladder"),
  c("TCGA-LGG", "Brain"),
  c("TCGA-BRCA", "Breast"),
  c("TCGA-CESC", "Cervix Uteri"),
  c("TCGA-CHOL", "Liver"), # ??
  c("TCGA-COAD", "Colon"),
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
  c("TCGA-UCEC", "Uterus")
  #c("TCGA-UVM", "") # ??
))
colnames(MATCHES) <- c("TCGA", "GTEX")

all_data <- load_data("~/Files/data/tcga/TCGA/") 
data <- clean_data(all_data)
rm(all_data)

gtex_meta <- read_tsv("./GTEX_phenotype")
gtex_data <- read_tsv("~/Files/data/tcga/gtex_gene_expected_count")
gtex_data_split <- split_gtex_data(gtex_data, gtex_meta)
rm(gtex_data)
gtex_data_split <- clean_data(gtex_data_split, "sample")

run_deseq(data, gtex_data_split, MATCHES, "~/Files/data/tcga/TCGA_dea_output/")
