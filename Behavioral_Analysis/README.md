# Behavioral Analysis

This repository provides scripts and resources used for behavioral analysis in [this study](https://www.biorxiv.org/content/10.1101/2024.12.20.629776v1).

---

## Environment

All analyses were performed using the following software versions:

- **Python**: 3.11.5  
- **Jupyter lab**: 3.6.3  

### Required Python Packages

```
numpy       # ver 1.24.3
pandas      # ver 2.0.3
matplotlib  # ver 3.7.2
seaborn     # ver 0.12.2
scipy       # ver 1.11.3
```

---

## DeepLabCut (DLC) Analysis

DeepLabCut (DLC) **v2.3.8** was used for pose estimation–based behavioral analysis.

- Official installation guide:  
  https://deeplabcut.github.io/DeepLabCut/

### Notebooks

- **`1_DLC_analysis_create_csv_v1.ipynb`**  
  Processes DeepLabCut output files and generates CSV summary files for downstream analysis.

- **`2_DLC_analysis_create_figure_v1.ipynb`**  
  Generates figures and plots using the CSV summaries produced by  
  `1_DLC_analysis_create_csv_v1.ipynb`.

### Required Python Scripts

- [`DLC_sqOF_v1.py`](./Scripts/DLC_sqOF_v1.py)
- [`DLC_cyOF_v1.py`](./Scripts/DLC_cyOF_v1.py)
- [`DLC_plot_v1.py`](./Scripts/DLC_plot_v1.py)

---

## Motion Sensor Analysis
### Notebook

A dedicated notebook is provided for **motion sensor–based behavioral analysis**:

- **`Motion_sensor_analysis_v1.ipynb`**

### Hardware

The motion sensor devices used in this study are part of the [**Harp** platform](https://www.cf-hw.org/harp/wear):

- **WEAR Basestation**
- **WEAR – Wireless sensor device**

---

## Notes

- Scripts are versioned (`_v1`) to support reproducibility and future updates.

---

## Citation

If you use this code, please cite the corresponding publication (details to be added).

---

## Contributions

DLC analysis scripts were developed by Akira Fushiki (Columbia University/Allen Institute). Motion sensor analysis scripts were primarily developed by Joaquim Alves da Silva (Champalimaud Foundation), with adjustments made by Akira Fushiki.
Feedback and contributions are welcome.

If you identify any errors, have suggestions for improvement, or wish to contribute to the codebase, please open an issue or submit a pull request via GitHub. All contributions that improve clarity, functionality, or reproducibility are appreciated.

