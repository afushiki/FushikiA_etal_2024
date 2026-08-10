# Sequencing Analysis

Scripts for the single-nucleus RNA-sequencing (snRNA-seq) analyses in
[Fushiki et al.](https://www.biorxiv.org/content/10.1101/2024.12.20.629776v2) — quality control,
ambient RNA removal, doublet detection, cell type annotation, differential abundance testing, and
downstream statistical and functional analyses.

← [Back to repository overview](../README.md)

---

## Contents

```
Sequencing_Analysis/
├── Scripts/
│   ├── Sequencing_analysis_v1.R
│   ├── SoupX_function_v1.R
│   ├── Doublet_function_v1.R
│   ├── MASC_analysis_v1.R
│   ├── Milo_analysis_v1.R
│   ├── Proportion_analysis_df_stat_v1.R
│   ├── Proportion_analysis_plot_v1.R
│   └── GO_analysis_v1.R
└── Annotation/
    └── mapmycells-output-af_1715006792868.csv
```

---

## Environment

| Software | Version |
| --- | --- |
| R | 4.3.3 |
| RStudio | 2023.12.1 |

<details>
<summary>R package versions</summary>

```r
Seurat                # v5.0.3
tidyverse             # v2.0.0
cowplot               # v1.1.3
DoubletFinder         # v2.0.4
SeuratWrappers        # v0.3.5
reticulate            # v1.35.0
patchwork             # v1.2.0
scCustomize           # v2.1.2
SingleCellExperiment  # v1.24.0
miloR                 # v1.99.9
lme4                  # v1.1-35.3
anndata               # v0.7.5.6
gprofiler2            # v0.2.3
gt                    # v0.10.1
gtsummary             # v1.7.2
SoupX                 # v1.6.2
DropletUtils          # v1.22.0
knitr                 # v1.46
RColorBrewer          # v1.1-3
ggpubr                # v0.6.0
```

</details>

---

## Data

Sequencing data are deposited in the Gene Expression Omnibus
([GSE285508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285508)).

Nuclei were isolated following
[this protocol](https://www.protocols.io/view/isolation-of-nuclei-from-adult-mouse-brain-tissue-14egn8b36g5d/v2).

---

## Pipeline

[`Sequencing_analysis_v1.R`](./Scripts/Sequencing_analysis_v1.R) is the entry point: it runs
preprocessing, quality control, clustering, and visualization, sourcing the modules below as needed.

**Preprocessing**

| Script | Purpose | Upstream method |
| --- | --- | --- |
| [`SoupX_function_v1.R`](./Scripts/SoupX_function_v1.R) | Estimates and removes ambient RNA contamination | [SoupX](https://github.com/constantAmateur/SoupX) |
| [`Doublet_function_v1.R`](./Scripts/Doublet_function_v1.R) | Identifies and removes doublets | [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder) |

**Differential abundance and composition**

| Script | Purpose | Upstream method |
| --- | --- | --- |
| [`MASC_analysis_v1.R`](./Scripts/MASC_analysis_v1.R) | Mixed-effects modeling of associations of single cells, testing cell population associations with disease status | [MASC](https://github.com/immunogenomics/masc) |
| [`Milo_analysis_v1.R`](./Scripts/Milo_analysis_v1.R) | Differential abundance testing on KNN graphs | [miloR](https://github.com/MarioniLab/miloR) |
| [`Proportion_analysis_df_stat_v1.R`](./Scripts/Proportion_analysis_df_stat_v1.R) | Builds summary data frames for proportional analysis of clusters across samples | — |
| [`Proportion_analysis_plot_v1.R`](./Scripts/Proportion_analysis_plot_v1.R) | Plots cluster proportions with statistical comparisons across conditions | — |

**Functional analysis**

| Script | Purpose | Upstream method |
| --- | --- | --- |
| [`GO_analysis_v1.R`](./Scripts/GO_analysis_v1.R) | Gene Ontology enrichment analysis, with tabular and graphical output | [gprofiler2](https://github.com/egonw/r-gprofiler2) |

---

## Cell type annotation

Cell types were annotated with **MapMyCells**, developed by the NIH BRAIN Initiative Cell Census
Network in collaboration with the Allen Institute
(<https://portal.brain-map.org/atlases-and-data/bkp/mapmycells>).

| Setting | Value |
| --- | --- |
| Reference taxonomy | 10x Whole Mouse Brain (CCN20230722) |
| Mapping algorithm | Hierarchical mapping |
| Output | [`mapmycells-output-af_1715006792868.csv`](./Annotation/mapmycells-output-af_1715006792868.csv) |

The output CSV is read back into the pipeline to label clusters.

---

## Notes

- Filenames are versioned (`_v1`) so that future revisions can be added without breaking references
  in the manuscript or in downstream analyses.
- External methods are credited via their original repositories; please cite them alongside this work.
- The pipeline is modular and can be adapted to other snRNA-seq datasets with minimal modification.

---

## Authors

Sequencing analysis scripts were written primarily by Akira Fushiki
(Columbia University / Allen Institute), with guidance from Zachary Lewis (Allen Institute) and
Archana Yadav (Columbia University).

Corrections, suggestions, and contributions are welcome — please
[open an issue or pull request](https://github.com/afushiki/FushikiA_etal_2024/issues).

## Citation

See the [repository overview](../README.md#citation) for the full citation and BibTeX entry.
