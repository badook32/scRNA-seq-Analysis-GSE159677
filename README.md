# GSE159677 scRNA-seq Analysis

This repository contains an R script for analyzing single-cell RNA-seq (scRNA-seq) data from the GEO dataset [GSE159677](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159677).  
The workflow is based on the Seurat package and includes steps from data loading and quality control to clustering, visualization, and cell cycle analysis.

> **Note**: This script was adapted and customized based on the excellent scRNA-seq workflow by Roman Haase:  
> [https://romanhaa.github.io/projects/scrnaseq_workflow/](https://romanhaa.github.io/projects/scrnaseq_workflow/)

---

## 📁 File Overview

| File | Description |
|------|-------------|
| `GSE159677_Analysis.R` | Main R script for processing and analyzing the dataset |
| `plots/` | Directory where quality control and analysis plots are saved |
| `data/` | Directory to save intermediate processed data |

---

## 🔬 Analysis Workflow

1. **Package Loading and Color Palette Setup**
2. **10X Genomics Data Import**  
   - Using custom `load10xData()` function
3. **Gene Matrix Cleaning and Sample Merging**
4. **Cell Quality Control**  
   - Mitochondrial transcript percentage  
   - Doublet detection with `scDblFinder`  
   - Filtering by nCount, nFeature, and percent.mt
5. **Gene Filtering**  
   - Retain genes expressed in ≥5 cells
6. **Seurat Object Construction and Normalization**  
   - Using `LogNormalize` and optional `SCTransform`
7. **Dimensionality Reduction & Clustering**  
   - PCA, clustering with `FindClusters`
8. **Cluster Stability Assessment**  
   - Bootstrap-based probability matrix  
   - Silhouette score analysis  
   - Cluster modularity and tree construction
9. **Sample-Cluster Composition Comparison**
10. **Cell Cycle Scoring**

---

## 🛠 Dependencies

The script uses the following R packages:

- `Seurat`
- `scran`, `scater`, `scDblFinder`, `SingleR`
- `SingleCellExperiment`, `bluster`, `cluster`, `cerebroApp`
- `ggplot2`, `patchwork`, `ggforce`, `viridis`, `ggridges`

Install all dependencies using:

```r
install.packages("tidyverse")
BiocManager::install(c("Seurat", "SingleCellExperiment", "scDblFinder", "SingleR", "bluster", "cerebroApp", "scran", "scater", "cluster"))


▶️ How to Run
Clone this repository or download the R script.

Ensure your working directory is set to the folder containing the GSE159677 dataset in 10X Genomics format.
setwd("~/Desktop/.../GSE159677")
source("GSE159677_Analysis.R")


📄 License
This repository is released under the MIT License.
However, please also review the license of the referenced tutorial:
https://romanhaa.github.io/projects/scrnaseq_workflow/


📚 References
https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE159677
https://romanhaa.github.io/projects/scrnaseq_workflow/
https://satijalab.org/seurat/
