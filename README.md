# Differential Gene Expression (DGE) Analysis with DESeq2: A Step-by-Step Workflow
A reproducible pipeline for Differential Gene Expression analysis using DESeq2, complete with step-by-step documentation, example data, and ready-to-use scripts in R.

## Purpose of This Repository
This repository provides an accessible, reproducible, and well-documented pipeline for Differential Gene Expression (DGE) analysis.
1. Installation and setup
2. Data preparation and preprocessing
3. Running differential expression analysis
4. Exporting normalized counts and results
5. Data transformation for visualization
6. Combining results for interpretation

## Installation & Setup
To install all the required R packages, execute the following commands in your R environment:
```r
if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
BiocManager::install(c("DESeq2", "apeglm"))
library(DESeq2)
library(apeglm)
```

## 1. Data Preparation
Prepare two essential input files:
- counts.csv: Raw counts of gene expression (genes in rows, samples in columns).
- sample_info.txt: Sample metadata including experimental conditions.

Load the data into R:
```r
count_data <- read.csv("counts.csv", header = TRUE, row.names = 1)
sample_info <- read.table("sample_info.txt", header = TRUE, sep = '\t')
```

## 2. Data Preprocessing
Replace zeros with NA and remove genes with incomplete data:
```r
count_data[count_data == 0] <- NA
filtered_counts <- count_data[complete.cases(count_data), ]
```

## 3. Running DESeq2 Analysis
Create a DESeq2 dataset object and perform the analysis:
```r
dds <- DESeqDataSetFromMatrix(countData = filtered_counts, 
                              colData = sample_info, 
                              design = ~ condition)

# Filter out lowly expressed genes
keep_genes <- rowSums(counts(dds)) >= 10
dds <- dds[keep_genes, ]

# Run DESeq2 analysis
dds <- DESeq(dds)
```

## 4. Results and Outputs
Export normalized count data for downstream analysis:
```r
normalized_counts <- counts(dds, normalized = TRUE)
write.csv(normalized_counts, "normalized_counts.csv")
```

Export Differential Expression Results
```r
results_dds <- results(dds, alpha = 0.05)
significant_results <- results_dds[order(results_dds$padj), ]
write.csv(significant_results, "DEG_results.csv")
```

## 5. Visualization Preparation
Transform data for visualization using regularized log transformation:
```r
rlog_transformed <- rlog(dds)
write.csv(assay(rlog_transformed), file = "rlog_transformed_counts.csv", quote = FALSE, row.names = TRUE)
```
6. Combining Results for Interpretation
Merge differential expression results with normalized counts for comprehensive analysis:
```r
final_data <- merge(as.data.frame(results_dds), 
                    as.data.frame(normalized_counts), 
                    by = "row.names", sort = FALSE)
names(final_data)[1] <- "Gene"
write.csv(final_data, file = "final_results.csv", quote = FALSE, row.names = FALSE)
```

## Output Summary
The analysis generates the following files:

File	Description
normalized_counts.csv	Normalized counts data for further analyses
DEG_results.csv	Sorted list of significantly differentially expressed genes
rlog_transformed_counts.csv	Regularized log-transformed counts for visualization
final_results.csv	Comprehensive combined results and normalized counts
![image](https://github.com/user-attachments/assets/ee06d7e3-2ef4-4728-8e9a-bf9ae3460e6c)

📬 Contact
For questions, suggestions, or issues, please reach out via:

Email: bioinfosourabh@gmail.com


