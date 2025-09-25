# B4galt1 Data Analysis
## B4galt1 single-cell analysis (LSK)

This repo contains code and example data to demo and reproduce the B4galt1 LSK single‑cell RNA‑seq analysis. Main workflow: `B4galt1_publication.Rmd` includes **transfer learning** (ML) for cell type classification. Optional scripts: `mk_velocyto.R` (RNA velocity) and `trajectory.R` (Monocle3).

## 1) System requirements
- OS: macOS (tested), Linux (expected compatible)
- CPU: Standard desktop/laptop; no GPU required
- RAM: 16+ GB (core); 32+ GB recommended for velocity/trajectory
- Disk: 20+ GB free

Tested software versions:
- R 4.2–4.3
- Key R packages: Seurat 4.3.0, sctransform 0.4.0, SeuratWrappers, SingleR, celldex, scater, tricycle, cowplot, patchwork, ggplot2, ggrepel, ggpubr, aplot, gridExtra, dplyr, tidyr, data.table, cluster, diptest, **xgboost** (for transfer learning ML), pivottabler, openxlsx, clusterProfiler, org.Mm.eg.db, msigdbr, velocyto.R (optional), monocle3 (optional).

No non‑standard hardware required.

## 2) Installation guide
1. Install R (4.3 recommended) and optionally RStudio.
2. Install packages:
```r
install.packages(c(
  "Seurat","SeuratWrappers","sctransform","SingleR","celldex","scater",
  "tricycle","cowplot","patchwork","ggplot2","ggrepel","ggpubr","aplot",
  "gridExtra","dplyr","tidyr","data.table","cluster","diptest","xgboost",
  "pivottabler","openxlsx","clusterProfiler","msigdbr"
))
if (!requireNamespace("BiocManager", quietly=TRUE)) install.packages("BiocManager")
BiocManager::install("org.Mm.eg.db")
# Optional
install.packages("velocyto.R")
BiocManager::install("monocle3")
```
Typical install time: 10–30 min on a normal desktop.

## 3) Demo
A minimal 10x‑style dataset is provided under `data/`.

Run main analysis:
```r
rmarkdown::render("B4galt1_publication.Rmd")
```
Outputs are written into an auto‑named folder (starts with the sample list and `mito_20_seurat4.3.0_sct0.4.0_ccGenes/`) containing PDFs and RDS files (e.g., `combined_SCT_finalClusters.Rds`). The analysis includes **transfer learning** using XGBoost for cell type classification (LTHSC, STHSC, MPP2, MPP3, MPP4).

Expected runtime: 60–120 min for the main Rmd on a normal desktop.

RNA velocity (optional): adjust paths inside `mk_velocyto.R` to your loom files and combined object, then:
```r
source("mk_velocyto.R")
```
Expected runtime: ~30–90 min.

Trajectory (optional): use the combined object from main script, then:
```r
source("trajectory.R")
```
Expected runtime: ~20–60 min.

## 4) Instructions for use (your data)
- Place your 10x outputs (barcodes/features/matrix) under a directory like `data/`, one subfolder per sample.
- In `B4galt1_publication.Rmd`, update:
  - `samples` and `phenotypes`
  - Path in `input_10x(...)`
  - `mito_cutoff` as desired
- Knit the Rmd to generate integrated objects and figures.
- Transfer learning data are under `data/` folder. change path to locate correct files. 

Notes:
- If adding/removing samples (e.g., removed `LSK_KO_2`), ensure any hardcoded lines referencing that sample are removed or replaced by loops.
- Update file paths in `mk_velocyto.R` and `trajectory.R` to your environment before running.

## (Optional) Reproduction instructions
- Knit `B4galt1_publication.Rmd` to regenerate all main figures and RDS outputs under the output folder.
- Then (optionally) run `mk_velocyto.R` and `trajectory.R` after path adjustments to reproduce velocity/trajectory outputs.
- Session info files (e.g., `sessionInfo.txt`, `finalsessionInfo.txt`) are produced to capture versions.

## Contact
Questions about running or reproducing results: contact Hoffmeister lab maintainers.
