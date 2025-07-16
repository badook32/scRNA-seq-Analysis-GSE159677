# GSE159677 scRNA-seq Analysis

This repository contains an R script for analyzing single-cell RNA-seq (scRNA-seq) data from the GEO dataset [GSE159677](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159677).  
The workflow is based on the Seurat package and includes steps from data loading and quality control to clustering, visualization, and cell cycle analysis.

> **Note**: This script was adapted and customized based on the excellent scRNA-seq workflow by Roman Haase:  
> [https://romanhaa.github.io/projects/scrnaseq_workflow/](https://romanhaa.github.io/projects/scrnaseq_workflow/)

---

## 🔬 Analysis Workflow

1. **Package Loading and Color Palette Setup**
2. **10X Genomics Data Import**  
   - Using custom `load10xData()` function
3. **Gene Matrix Cleaning and Sample Merging**
4. **Cell Quality Control**  
   - Mitochondrial transcript percentage  
   - Doublet detection with `scDblFinder`  
   - Filtering by `nCount`, `nFeature`, and `percent.mt`
5. **Gene Filtering**  
   - Retain genes expressed in ≥5 cells
6. **Seurat Object Construction and Normalization**  
   - Using `LogNormalize` and optionally `SCTransform`
7. **Dimensionality Reduction & Clustering**  
   - PCA, UMAP, and `FindClusters`
8. **Cluster Stability Assessment**  
   - Bootstrap-based probability matrix  
   - Silhouette score analysis  
   - Modularity and cluster tree visualization
9. **Sample-Cluster Composition Comparison**
10. **Cell Cycle Scoring**

---

## 🛠 Dependencies

This script requires the following R packages:

- `Seurat`
- `scran`, `scater`, `scDblFinder`, `SingleR`
- `SingleCellExperiment`, `bluster`, `cluster`, `cerebroApp`
- `ggplot2`, `patchwork`, `ggforce`, `viridis`, `ggridges`

Install all dependencies using:

```r
# CRAN packages
install.packages(c("tidyverse", "ggforce", "ggridges", "patchwork", "viridis", "cluster"))

# Bioconductor manager
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

# Bioconductor packages
BiocManager::install(c(
  "Seurat",
  "SingleCellExperiment",
  "scDblFinder",
  "SingleR",
  "bluster",
  "cerebroApp",
  "scran",
  "scater"
))
```

---
## ▶️ How to Run

Clone this repository or download the R script file.

Ensure your working directory is set to the folder containing the GSE159677 dataset (in 10X Genomics format), and then run the script:

```r
setwd("~/Desktop/.../GSE159677")
source("GSE159677_Analysis.R")
```

---
## 📚 References
GSE159677 - NCBI GEO
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159677

Roman Haase's scRNA-seq workflow
https://romanhaa.github.io/projects/scrnaseq_workflow/

Seurat documentation
https://satijalab.org/seurat/
