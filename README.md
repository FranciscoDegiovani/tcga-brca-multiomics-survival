# Multi-Omics Survival Analysis — TCGA-BRCA (MOFA2)

Integration of RNA-seq and DNA methylation data from breast cancer patients (TCGA-BRCA) using MOFA2 (Multi-Omics Factor Analysis) to identify latent factors, followed by Cox regression and Kaplan-Meier analysis to test their association with patient survival. Built as part of my bioinformatics/genomics portfolio, applying multi-omics integration methodology — the same type of approach used in my master's research — to a fully public cohort and disease area unrelated to my thesis (breast cancer here, vs. endometrial cancer in my actual research).

## Data Source

- **Cohort:** TCGA-BRCA (The Cancer Genome Atlas — Breast Invasive Carcinoma), accessed via the NCI Genomic Data Commons (GDC)
- **Data types:** RNA-seq gene expression (STAR - Counts), DNA methylation (Illumina 450K array), clinical/survival data
- **Access tool:** TCGAbiolinks (R/Bioconductor)
- **Citation:** Colaprico A, Silva TC, Olsen C, et al. (2016). TCGAbiolinks: an R/Bioconductor package for integrative analysis of TCGA data. *Nucleic Acids Research*, 44(8), e71.
- **License:** TCGA data is open-access (tumor/normal, non-controlled tier), publicly available for research use

## Tools

- R / Bioconductor (TCGAbiolinks, SummarizedExperiment, survival, survminer)
- Python (mofapy2) via Google Colab — used for MOFA2 training due to a `basilisk`/Python environment installation issue on Windows
- MOFA2 (Multi-Omics Factor Analysis) — Argelaguet R, et al. (2018). Multi-Omics Factor Analysis—a framework for unsupervised integration of multi-omics data sets. *Molecular Systems Biology*, 14(6), e8124.

## Workflow

1. Query and download RNA-seq and DNA methylation data for matched tumor/normal samples (TCGAbiolinks/GDC)
2. Manually assemble the methylation beta-value matrix (bypassing `sesame`, unavailable for the Bioconductor version used)
3. Match samples present in both omics layers by patient barcode (866 matched samples)
4. Normalize RNA-seq (log2 CPM) and select the top 5,000 most variable features per omic layer
5. Train a MOFA2 model (10 factors) in Python, using the same preprocessed data exported from R
6. Load the trained model back into R and extract variance explained per view/factor
7. Query clinical/survival data and merge with MOFA2 factors by patient ID (785 matched patients)
8. Test each factor's association with survival via Cox proportional hazards regression
9. Visualize the two strongest factors with Kaplan-Meier curves

## Results

- **Factor 1** explains the largest share of variance, especially in methylation (~30%)
- **Cox regression** identified two factors significantly associated with survival: **Factor4** (p = 0.0004, HR = 0.87) and **Factor5** (p = 0.0034, HR = 0.89) — both act as protective factors (higher values associated with lower risk of death); **Factor10** was also significant (p = 0.037)
- **Kaplan-Meier (Factor4):** clear separation between high/low groups, log-rank p = 0.00016
- **Kaplan-Meier (Factor5):** log-rank p = 0.081 (not significant at the median split), despite a highly significant continuous Cox result (p = 0.0034) — a good illustration of how median-splitting a continuous variable loses statistical power compared to modeling it continuously

## Structure

scripts/
01_setup_and_packages.R (accumulates all preprocessing and analysis steps)
notebooks/
mofa2_training.ipynb
results/
figures/
km_factor4.png
km_factor5.png
tables/
tcga_brca_rna_query_results.csv
tcga_brca_meth_query_results.csv
cox_regression_by_factor.csv


Note: large intermediate files (trained MOFA2 model, full filtered expression/methylation matrices) are excluded from version control (see `.gitignore`) and kept local only — standard practice for reproducibility without repository bloat.

## Skills Demonstrated

- Multi-omics data acquisition and integration (GDC/TCGAbiolinks)
- Cross-platform workflow (R for preprocessing/statistics, Python for model training) with careful sample ID matching across the pipeline
- Unsupervised multi-omics factor analysis (MOFA2)
- Survival analysis: Cox proportional hazards regression and Kaplan-Meier estimation
- Debugging real-world data integration issues (ID format mismatches, environment/dependency errors)
- Git version control discipline (excluding large data artifacts from the repository)

## Status

✅ Complete
