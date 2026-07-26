# 01_setup_and_packages.R
# Load required packages for the TCGA-BRCA multi-omics survival project

library(TCGAbiolinks)
library(MOFA2)
library(survival)
library(survminer)

# 02_query_rna_data.R
# Query TCGA-BRCA RNA-seq (STAR - Counts) data available on GDC.
# This step only inspects what's available — no download yet.

library(TCGAbiolinks)

query_rna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

rna_results <- getResults(query_rna)
dim(rna_results)
head(rna_results[, c("cases", "sample_type")])
table(rna_results$sample_type)

# 02_query_rna_data.R
# Query TCGA-BRCA RNA-seq (STAR - Counts) data available on GDC.
# This step only inspects what's available — no download yet.

library(TCGAbiolinks)

query_rna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

rna_results <- getResults(query_rna)
dim(rna_results)
head(rna_results[, c("cases", "sample_type")])
table(rna_results$sample_type)

# Save the query results table for reference/documentation
dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)
write.csv(rna_results, "results/tables/tcga_brca_rna_query_results.csv", row.names = FALSE)

# 03_download_rna_data.R
# Download RNA-seq counts for tumor and normal samples only
# (excluding the 7 metastatic samples, too few for robust separate analysis)

library(TCGAbiolinks)

query_rna <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

GDCdownload(query_rna, directory = "data/raw")

# Confirm the download completed and check the local file structure
list.files("data/raw", recursive = FALSE)
length(list.files("data/raw", recursive = TRUE, pattern = "\\.tsv$"))


# 04_prepare_rna_matrix.R
# Consolidate the downloaded per-sample files into a single
# SummarizedExperiment object (genes x samples expression matrix
# with clinical metadata attached).

library(TCGAbiolinks)

rna_data <- GDCprepare(query = query_rna, directory = "data/raw")

dim(rna_data)
class(rna_data)

saveRDS(rna_data, "data/processed/tcga_brca_rna_raw.rds")


