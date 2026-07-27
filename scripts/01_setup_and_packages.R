# 01_setup_and_packages.R
# Load required packages for the TCGA-BRCA multi-omics survival project

library(TCGAbiolinks)
library(MOFA2)
library(survival)
library(survminer)
library(SummarizedExperiment)

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


# 08_match_samples_and_prepare_mofa_input.R
# Match samples present in both RNA-seq and methylation data,
# then prepare the two matrices as input views for MOFA2.

library(TCGAbiolinks)

rna_data <- readRDS("data/processed/tcga_brca_rna_raw.rds")
meth_matrix <- readRDS("data/processed/tcga_brca_methylation_matrix.rds")

# TCGA barcodes: RNA-seq colnames are full barcodes, methylation colnames are file UUIDs
# Need to map methylation UUIDs back to sample barcodes using the query results table
meth_results <- read.csv("results/tables/tcga_brca_rna_query_results.csv")  # placeholder check

# Get sample barcode for each RNA-seq sample (first 15 characters = patient + sample type)
rna_barcodes <- substr(colnames(rna_data), 1, 15)

# Check overlap
length(rna_barcodes)
length(colnames(meth_matrix))


# Save methylation query results (was missing since step 05)
meth_query_results <- getResults(query_meth)
write.csv(meth_query_results, "results/tables/tcga_brca_meth_query_results.csv", row.names = FALSE)

# Check which columns map UUID (id) to sample barcode (cases)
head(meth_query_results[, c("id", "cases")])


# Map methylation UUIDs (column names) to sample barcodes
meth_id_to_barcode <- setNames(meth_query_results$cases, meth_query_results$id)

meth_barcodes_full <- meth_id_to_barcode[colnames(meth_matrix)]
sum(is.na(meth_barcodes_full))  # should be 0 — check every UUID found a match

# Truncate to first 15 characters, same format used for RNA-seq matching
meth_barcodes <- substr(meth_barcodes_full, 1, 15)

# Find samples present in BOTH omics
common_barcodes <- intersect(rna_barcodes, meth_barcodes)
length(common_barcodes)


rna_counts <- assay(rna_data)

# Subset RNA-seq matrix to common samples (handle possible duplicate barcodes)
rna_counts <- assay(rna_data)
colnames(rna_counts) <- rna_barcodes

# In case of duplicate barcodes in RNA-seq, keep only the first occurrence
rna_counts <- rna_counts[, !duplicated(colnames(rna_counts))]

rna_common <- rna_counts[, common_barcodes]

# Subset methylation matrix the same way
colnames(meth_matrix) <- meth_barcodes
meth_matrix <- meth_matrix[, !duplicated(colnames(meth_matrix))]

meth_common <- meth_matrix[, common_barcodes]

# Confirm both matrices are now aligned
dim(rna_common)
dim(meth_common)
identical(colnames(rna_common), colnames(meth_common))

rna_counts <- assay(rna_data)
colnames(rna_counts) <- rna_barcodes
rna_counts <- rna_counts[, !duplicated(colnames(rna_counts))]
rna_common <- rna_counts[, common_barcodes]

colnames(meth_matrix) <- meth_barcodes
meth_matrix <- meth_matrix[, !duplicated(colnames(meth_matrix))]
meth_common <- meth_matrix[, common_barcodes]

dim(rna_common)
dim(meth_common)
identical(colnames(rna_common), colnames(meth_common))


# 09_normalize_and_select_features.R
# Normalize RNA-seq counts (log2 CPM) and select the top variable
# features from each omic layer before MOFA2.

# RNA-seq: convert raw counts to log2(CPM + 1)
lib_sizes <- colSums(rna_common)
rna_cpm <- t(t(rna_common) / lib_sizes) * 1e6
rna_log <- log2(rna_cpm + 1)

# Select top 5000 most variable genes
rna_var <- apply(rna_log, 1, var)
top_rna_genes <- names(sort(rna_var, decreasing = TRUE))[1:5000]
rna_filtered <- rna_log[top_rna_genes, ]

# Methylation: already beta values (0-1 scale), no need to normalize
# Select top 5000 most variable CpG probes
# (using a subsample of rows first to estimate variance faster, then filtering)
meth_var <- apply(meth_common, 1, var)
top_meth_probes <- names(sort(meth_var, decreasing = TRUE))[1:5000]
meth_filtered <- meth_common[top_meth_probes, ]

dim(rna_filtered)
dim(meth_filtered)


# 10_run_mofa2.R
# Build and train the MOFA2 multi-omics model

library(MOFA2)

mofa_data <- list(
  RNA = as.matrix(rna_filtered),
  Methylation = as.matrix(meth_filtered)
)

mofa_object <- create_mofa(mofa_data)

# Inspect the data overview before training
plot_data_overview(mofa_object)

# Set training options
data_opts <- get_default_data_options(mofa_object)
model_opts <- get_default_model_options(mofa_object)
train_opts <- get_default_training_options(mofa_object)

# Set number of factors to learn (10 is a common starting point)
model_opts$num_factors <- 10

# Reduce verbosity but keep progress visible; cap iterations for reasonable runtime
train_opts$convergence_mode <- "fast"
train_opts$seed <- 42

mofa_object <- prepare_mofa(
  object = mofa_object,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

mofa_train

# 10b_export_matrices_for_python.R
# Export filtered matrices as CSV for training MOFA2 in Python (Colab),
# due to basilisk/pyenv installation issues on Windows.

dir.create("results/tables", showWarnings = FALSE, recursive = TRUE)

# Transpose so rows = samples, columns = features (standard format for mofapy2)
write.csv(t(rna_filtered), "results/tables/rna_filtered_samples_by_features.csv")
write.csv(t(meth_filtered), "results/tables/meth_filtered_samples_by_features.csv")


mofa_trained <- load_model("results/mofa_model.hdf5")
plot_variance_explained(mofa_trained, x = "view", y = "factor")


# 11_get_clinical_data.R
# Query clinical/survival data for the same TCGA-BRCA patients

library(TCGAbiolinks)

clinical_data <- GDCquery_clinic(project = "TCGA-BRCA", type = "clinical")

dim(clinical_data)
head(clinical_data[, c("submitter_id", "vital_status", "days_to_death", "days_to_last_follow_up")])



# 12_prepare_survival_data.R
# Build the time-to-event variable and match clinical data to
# the MOFA2 factors by patient ID.

library(dplyr)

# Combine days_to_death (if dead) or days_to_last_follow_up (if alive)
# into a single "time" variable, and convert vital status to numeric event
clinical_surv <- clinical_data %>%
  mutate(
    time = ifelse(vital_status == "Dead", days_to_death, days_to_last_follow_up),
    event = ifelse(vital_status == "Dead", 1, 0)
  ) %>%
  filter(!is.na(time)) %>%
  select(submitter_id, time, event)

dim(clinical_surv)
head(clinical_surv)

# Extract MOFA2 factors
factors <- get_factors(mofa_trained, factors = "all")[[1]]

# Match factors' sample barcodes (15-char) to clinical patient IDs (12-char)
patient_ids <- substr(rownames(factors), 1, 12)
rownames(factors) <- patient_ids

# Merge factors with survival data
factors_df <- as.data.frame(factors)
factors_df$submitter_id <- rownames(factors_df)

survival_merged <- inner_join(factors_df, clinical_surv, by = "submitter_id")

dim(survival_merged)
head(survival_merged)


# Diagnostic: compare ID formats before and after truncation
head(rownames(get_factors(mofa_trained, factors = "all")[[1]]), 3)
head(patient_ids, 3)
head(clinical_surv$submitter_id, 3)

# Check for exact string match
intersect(patient_ids, clinical_surv$submitter_id) |> length()


# Compare exact character properties, not just visual appearance
str(factors_df$submitter_id)
str(clinical_surv$submitter_id)

nchar(factors_df$submitter_id[1])
nchar(clinical_surv$submitter_id[1])

# Check for invisible characters (e.g. trailing spaces, encoding artifacts)
factors_df$submitter_id[1] == clinical_surv$submitter_id[1]

# Fix: replace dots back to hyphens in submitter_id
factors_df$submitter_id <- gsub("\\.", "-", factors_df$submitter_id)

head(factors_df$submitter_id, 3)

# Re-check overlap and re-merge
intersect(factors_df$submitter_id, clinical_surv$submitter_id) |> length()

survival_merged <- inner_join(factors_df, clinical_surv, by = "submitter_id")
dim(survival_merged)
head(survival_merged)


# 13_survival_analysis.R
# Test whether each MOFA2 factor is associated with patient survival
# using Cox proportional hazards regression.

library(survival)
library(survminer)

# Cox regression for each factor individually
cox_results <- lapply(paste0("Factor", 1:10), function(f) {
  formula <- as.formula(paste("Surv(time, event) ~", f))
  model <- coxph(formula, data = survival_merged)
  summary(model)$coefficients
})

names(cox_results) <- paste0("Factor", 1:10)

# Combine into a single results table
cox_summary <- do.call(rbind, lapply(names(cox_results), function(f) {
  data.frame(Factor = f, cox_results[[f]])
}))

cox_summary <- cox_summary[order(cox_summary$Pr...z..), ]
print(cox_summary)

write.csv(cox_summary, "results/tables/cox_regression_by_factor.csv", row.names = FALSE)


# 14_kaplan_meier_factor4.R
# Visualize the strongest factor (Factor4) by splitting patients
# into high/low groups at the median.

survival_merged$Factor4_group <- ifelse(
  survival_merged$Factor4 > median(survival_merged$Factor4), "High", "Low"
)

km_fit <- survfit(Surv(time, event) ~ Factor4_group, data = survival_merged)

km_plot <- ggsurvplot(
  km_fit,
  data = survival_merged,
  pval = TRUE,
  risk.table = TRUE,
  title = "Survival by Factor4 (MOFA2) — High vs. Low",
  xlab = "Days",
  legend.title = "Factor4 group",
  palette = c("firebrick", "steelblue")
)

png("results/figures/km_factor4.png", width = 900, height = 700)
print(km_plot)
dev.off()

km_plot


# 15_kaplan_meier_factor5.R
# Visualize the second-strongest factor (Factor5), same approach as Factor4.

survival_merged$Factor5_group <- ifelse(
  survival_merged$Factor5 > median(survival_merged$Factor5), "High", "Low"
)

km_fit_f5 <- survfit(Surv(time, event) ~ Factor5_group, data = survival_merged)

km_plot_f5 <- ggsurvplot(
  km_fit_f5,
  data = survival_merged,
  pval = TRUE,
  risk.table = TRUE,
  title = "Survival by Factor5 (MOFA2) — High vs. Low",
  xlab = "Days",
  legend.title = "Factor5 group",
  palette = c("firebrick", "steelblue")
)

png("results/figures/km_factor5.png", width = 900, height = 700)
print(km_plot_f5)
dev.off()

km_plot_f5

