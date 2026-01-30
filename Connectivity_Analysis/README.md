
# Connectivity Analysis

This repository provides scripts and resources used for connectivity analysis in [this study](https://www.biorxiv.org/content/10.1101/2024.12.20.629776v1).

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
```

---

### Data
- The [ABA_CCF_updated folder](./ABA_CCF_updated) contains revised annotations for the striatal subdivisions.

---

## Connectivity Analysis

### Notebooks

- [**`MitoPark_cell_analysis_v1.ipynb`**](./MitoPark_cell_analysis_v1.ipynb)  
  Processes MitoPark cell composition data and generate a CSV file and accompanying figure.

- [**`MitoPark_projection_analysis_v1.ipynb`**](./MitoPark_projection_analysis_v1.ipynb)  
  Processes MitoPark projection data and generate a CSV file and accompanying figure.

- [**`Connectivity_analysis_v1.ipynb`**](./Connectivity_analysis_v1.ipynb)  
  Processes [BrainJ](https://github.com/lahammond/BrainJ) outputs for rabies and projection datasets, producing summary CSV files and corresponding figures.


### Required Python Scripts

- [`brainJ_output_analysis_v1.py`](./brainJ_output_analysisv1.py)

---

## Notes

- Scripts are versioned (`_v1`) to support reproducibility and future updates.

---

## Citation

If you use this code, please cite the corresponding publication (details to be added).

---

## Contributions

Connectivity analysis scripts were developed by Akira Fushiki (Columbia University/Allen Institute). Feedback and contributions are welcome.

If you identify any errors, have suggestions for improvement, or wish to contribute to the codebase, please open an issue or submit a pull request via GitHub. All contributions that improve clarity, functionality, or reproducibility are appreciated.



