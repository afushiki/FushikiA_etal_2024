# Sequencing Analysis

This repository provides scripts and resources used for single-cell RNA-sequencing (scRNA-seq) analysis in [this study](https://www.biorxiv.org/content/10.1101/2024.12.20.629776v1), including quality control, doublet detection, ambient RNA removal, differential abundance testing, and downstream statistical and functional analyses.

---

## Environment

All analyses were performed using the following software versions:

- **R**: 4.3.3  
- **RStudio**: 2023.12.1  

### Required R Packages

```
library(Seurat)               # v5.0.3
library(tidyverse)            # v2.0.0
library(cowplot)              # v1.1.3
library(DoubletFinder)        # v2.0.4
library(SeuratWrappers)       # v0.3.5
library(reticulate)           # v1.35.0
library(patchwork)            # v1.2.0
library(scCustomize)          # v2.1.2
library(SingleCellExperiment) # v1.24.0
library(miloR)                # v1.99.9
library(lme4)                 # v1.1-35.3
library(anndata)              # v0.7.5.6
library(gprofiler2)           # v0.2.3
library(gt)                   # v0.10.1
library(gtsummary)            # v1.7.2
library(SoupX)                # v1.6.2
library(DropletUtils)         # v1.22.0
library(knitr)                # v1.46
library(RColorBrewer)         # v1.1-3
library(ggpubr)               # v0.6.0
```

---

## Data Availability

Single-cell RNA-sequencing data is deposited in the Gene Expression Omnibus website ([GSE285508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285508)). 

---

## Scripts Overview

The analysis pipeline is organized into modular R scripts:

1. **`Sequencing_analysis_v1.R`**  
   Primary script for preprocessing, quality control, clustering, and visualization of scRNA-seq data.

2. [**`SoupX_function_v1.R`**](./Scripts/SoupX_function_v1.R)  
   Estimates and removes ambient RNA contamination using *SoupX*.  
   Original implementation: https://github.com/constantAmateur/SoupX

3. [**`Doublet_function_v1.R`**](./Scripts/Doublet_function_v1.R)  
   Identifies and removes doublets using *DoubletFinder*.  
   Original implementation: https://github.com/chris-mcginnis-ucsf/DoubletFinder

4. **`MASC_analysis_v1.R`**  
   Performs MASC (Mixed-effects modeling of Associations of Single Cells) to test associations between cell populations and disease status at the single-cell level.  
   Original implementation: https://github.com/immunogenomics/masc

5. **`Milo_analysis_v1.R`**  
   Conducts differential abundance analysis using *miloR* on KNN graphs derived from scRNA-seq data.  
   Original implementation: https://github.com/MarioniLab/miloR

6. **`Proportion_analysis_df_stat_v1.R`**  
   Generates summary data frames for proportional analysis of cell clusters across samples.

7. **`Proportion_analysis_plot_v1.R`**  
   Produces proportional plots with statistical comparisons of cluster distributions across experimental conditions.

8. [**`GO_analysis_v1.R`**](./Scripts/GO_analysis_v1.R)  
   Performs Gene Ontology (GO) enrichment analysis using *g:Profiler* and generates tabular and graphical outputs.  
   Original implementation: https://github.com/egonw/r-gprofiler2

---

## Cell Type Annotation

Cell type annotation was performed using **MapMyCells**, an interactive online tool developed by the NIH BRAIN Initiative Cell Census Network in collaboration with the Allen Institute.

- **Tool**: https://portal.brain-map.org/atlases-and-data/bkp/mapmycells  
- **Reference Taxonomy**: 10x Whole Mouse Brain (CCN20230722)  
- **Mapping Algorithm**: Hierarchical Mapping  
- **Output**: Annotated CSV file imported into the analysis pipeline [available here](./Annotation/mapmycells-output-af_1715006792868.csv)

---

## Notes

- Scripts are versioned (`_v1`) to support reproducibility and future updates.
- External tools and methods are credited via their original GitHub repositories.
- The pipeline is modular and can be adapted to other scRNA-seq datasets with minimal modification.

---

## Citation

If you use this code, please cite the corresponding publication (details to be added).

---

## Contributions

Sequence analysis scripts were primarily developed by Akira Fushiki (Columbia University/Allen Institute), with guidance from Zack Lewis (Allen Institute) and Archana Yadav (Columbia University).
Feedback and contributions are welcome.

If you identify any errors, have suggestions for improvement, or wish to contribute to the codebase, please open an issue or submit a pull request via GitHub. All contributions that improve clarity, functionality, or reproducibility are appreciated.

