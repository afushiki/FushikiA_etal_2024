# Connectivity Analysis

Scripts for the anatomical and connectivity analyses in
[Fushiki et al.](https://www.biorxiv.org/content/10.1101/2024.12.20.629776v2) —
quantification of dopaminergic cell loss in MitoPark mice, and of the axonal projections and
monosynaptic inputs of *Anxa1*+ neurons.

← [Back to repository overview](../README.md)

---

## Contents

```
Connectivity_Analysis/
├── ABA_CCF_updated/                        # Revised striatal subdivision annotations
├── MitoPark_cell_analysis_v1.ipynb
├── MitoPark_projection_analysis_v1.ipynb
├── Connectivity_analysis_v1.ipynb
└── brainJ_output_analysis_v1.py
```

---

## Environment

| Software | Version |
| --- | --- |
| Python | 3.11.5 |
| JupyterLab | 3.6.3 |

```
numpy==1.24.3
pandas==2.0.3
matplotlib==3.7.2
seaborn==0.12.2
```

---

## Data

[`ABA_CCF_updated/`](./ABA_CCF_updated) contains revised annotations for the striatal subdivisions
of the Allen Mouse Brain Common Coordinate Framework, used to assign projections to striatal
subregions.

Whole-brain imaging datasets are deposited in the
[Brain Image Library](https://doi.brainimagelibrary.org/doi/10.35077/g.1209).

---

## Notebooks

Each notebook is self-contained and produces a summary CSV file together with the corresponding
figure panels.

| Notebook | Purpose |
| --- | --- |
| [`MitoPark_cell_analysis_v1.ipynb`](./MitoPark_cell_analysis_v1.ipynb) | Quantifies cell composition in MitoPark mice across timepoints |
| [`MitoPark_projection_analysis_v1.ipynb`](./MitoPark_projection_analysis_v1.ipynb) | Quantifies projection density in MitoPark mice |
| [`Connectivity_analysis_v1.ipynb`](./Connectivity_analysis_v1.ipynb) | Processes [BrainJ](https://github.com/lahammond/BrainJ) outputs for the rabies-tracing and anterograde projection datasets |

Supporting module, imported by the notebooks:

| Module | Purpose |
| --- | --- |
| [`brainJ_output_analysis_v1.py`](./brainJ_output_analysis_v1.py) | Parsing and quantification of BrainJ registration output |

Associated protocols:
[Costa Lab Nanoject virus injections](https://www.protocols.io/view/costa-lab-nanoject-virus-injections-cu-kxygxy6j4l8j/v1) ·
[Histology](https://www.protocols.io/view/histology-fushiki-et-al-2024-5qpvo9o49v4o/v1)

---

## Notes

Filenames are versioned (`_v1`) so that future revisions can be added without breaking references
in the manuscript or in downstream analyses.

---

## Authors

Connectivity analysis scripts were written by Akira Fushiki
(Columbia University / Allen Institute).

Corrections, suggestions, and contributions are welcome — please
[open an issue or pull request](https://github.com/afushiki/FushikiA_etal_2024/issues).

## Citation

See the [repository overview](../README.md#citation) for the full citation and BibTeX entry.
