# Behavioral Analysis

Scripts for the behavioral analyses in
[Fushiki et al.](https://www.biorxiv.org/content/10.1101/2024.12.20.629776v2) —
markerless pose estimation of open field and cylinder tests, and wireless motion-sensor recordings.

← [Back to repository overview](../README.md)

---

## Contents

```
Behavioral_Analysis/
├── DLC_Analysis/
│   ├── 1_DLC_analysis_create_csv_v1.ipynb
│   ├── 2_DLC_analysis_create_figure_v1.ipynb
│   ├── DLC_sqOF_v1.py
│   ├── DLC_cyOF_v1.py
│   └── DLC_plot_v1.py
└── Motion_Sensor_Analysis/
    ├── Motion_sensor_analysis_v1.ipynb
    └── motion_analysis_functions.py
```

---

## Environment

| Software | Version |
| --- | --- |
| Python | 3.11.5 |
| JupyterLab | 3.6.3 |
| DeepLabCut | 2.3.8 |

```
numpy==1.24.3
pandas==2.0.3
matplotlib==3.7.2
seaborn==0.12.2
scipy==1.11.3
```

DeepLabCut installation is documented at <https://deeplabcut.github.io/DeepLabCut/>.

---

## Data

Raw behavioral datasets and experimental video recordings are deposited in BioStudies
([S-BSST2649](https://www.ebi.ac.uk/biostudies/studies/S-BSST2649)). Download these before running
the notebooks.

---

## Pose estimation analysis (`DLC_Analysis/`)

Run the notebooks in numerical order; the second consumes the output of the first.

| Step | Notebook | Purpose |
| --- | --- | --- |
| 1 | [`1_DLC_analysis_create_csv_v1.ipynb`](./DLC_Analysis/1_DLC_analysis_create_csv_v1.ipynb) | Parses DeepLabCut output and writes CSV summary files for downstream analysis |
| 2 | [`2_DLC_analysis_create_figure_v1.ipynb`](./DLC_Analysis/2_DLC_analysis_create_figure_v1.ipynb) | Generates figures and plots from those CSV summaries |

Supporting modules, imported by the notebooks:

| Module | Purpose |
| --- | --- |
| [`DLC_sqOF_v1.py`](./DLC_Analysis/DLC_sqOF_v1.py) | Square open field analysis |
| [`DLC_cyOF_v1.py`](./DLC_Analysis/DLC_cyOF_v1.py) | Cylinder test analysis |
| [`DLC_plot_v1.py`](./DLC_Analysis/DLC_plot_v1.py) | Shared plotting functions |

Associated protocol:
[Open field & cylinder test behavior assays](https://www.protocols.io/view/open-field-amp-cylinder-test-behavior-assays-eq2ly6ypegx9/v1)

---

## Motion sensor analysis (`Motion_Sensor_Analysis/`)

| File | Purpose |
| --- | --- |
| [`Motion_sensor_analysis_v1.ipynb`](./Motion_Sensor_Analysis/Motion_sensor_analysis_v1.ipynb) | Analysis of wireless motion-sensor recordings |
| [`motion_analysis_functions.py`](./Motion_Sensor_Analysis/motion_analysis_functions.py) | Supporting functions imported by the notebook |

**Hardware.** Recordings were made with the [Harp](https://www.cf-hw.org/harp/wear) platform:
the WEAR basestation and WEAR wireless sensor device.

---

## Notes

Filenames are versioned (`_v1`) so that future revisions can be added without breaking references
in the manuscript or in downstream analyses.

---

## Authors

DeepLabCut analysis scripts were written by Akira Fushiki (Columbia University / Allen Institute).
Motion-sensor analysis scripts were written primarily by Joaquim Alves da Silva
(Champalimaud Foundation), with modifications by Akira Fushiki.

Corrections, suggestions, and contributions are welcome — please
[open an issue or pull request](https://github.com/afushiki/FushikiA_etal_2024/issues).

## Citation

See the [repository overview](../README.md#citation) for the full citation and BibTeX entry.
