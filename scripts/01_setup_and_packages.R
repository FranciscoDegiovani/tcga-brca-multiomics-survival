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



# 05_query_methylation_data.R
# Query TCGA-BRCA DNA methylation (450K array) data, matched to the
# same tumor/normal samples used for RNA-seq.

library(TCGAbiolinks)

query_meth <- GDCquery(
  project = "TCGA-BRCA",
  data.category = "DNA Methylation",
  platform = "Illumina Human Methylation 450",
  data.type = "Methylation Beta Value",
  sample.type = c("Primary Tumor", "Solid Tissue Normal")
)

meth_results <- getResults(query_meth)
dim(meth_results)
table(meth_results$sample_type)


# 06_download_methylation_data.R
# Estimate size before downloading (methylation arrays are large)

library(TCGAbiolinks)

# GDCquery already ran; just inspect estimated size via query metadata
sum(meth_results$file_size) / 1e9  # approximate size in GB

# Download DNA methylation beta values for the 890 matched tumor/normal samples

library(TCGAbiolinks)

GDCdownload(query_meth, directory = "data/raw")



# 07_prepare_methylation_matrix.R
# Manually assemble the methylation beta-value matrix from the
# downloaded per-sample files (bypassing GDCprepare, which requires
# the sesame/sesameData packages not yet available for this
# Bioconductor version).

meth_files <- list.files("data/raw", pattern = "methylation_array.*betas.*\\.txt$",
                         recursive = TRUE, full.names = TRUE)
length(meth_files)
head(meth_files, 3)

# Peek at the raw content of a single methylation file before
# writing the full assembly script
readLines(meth_files[1], n = 5)


# 07_prepare_methylation_matrix.R
# Manually assemble the methylation beta-value matrix (CpG probes x samples)
# from the downloaded per-sample files, bypassing GDCprepare/sesame.

library(data.table)

meth_files <- list.files("data/raw", pattern = "methylation_array.*betas.*\\.txt$",
                         recursive = TRUE, full.names = TRUE)

file_uuids <- basename(dirname(meth_files))

read_one_sample <- function(path) {
  dt <- fread(path, header = FALSE, col.names = c("probe_id", "beta"))
  dt$beta
}

message("Reading first file to get probe order...")
first_dt <- fread(meth_files[1], header = FALSE, col.names = c("probe_id", "beta"))
probe_ids <- first_dt$probe_id

message("Reading remaining ", length(meth_files), " files...")
beta_list <- lapply(seq_along(meth_files), function(i) {
  if (i %% 100 == 0) message("  ...", i, " / ", length(meth_files))
  read_one_sample(meth_files[i])
})

beta_matrix <- do.call(cbind, beta_list)
rownames(beta_matrix) <- probe_ids
colnames(beta_matrix) <- file_uuids

message("Final matrix dimensions: ", nrow(beta_matrix), " x ", ncol(beta_matrix))

saveRDS(beta_matrix, "data/processed/tcga_brca_methylation_matrix.rds")

message("Done. Saved to data/processed/tcga_brca_methylation_matrix.rds")