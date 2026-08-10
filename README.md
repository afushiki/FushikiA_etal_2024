# *Anxa1*+ Dopaminergic Neurons and Early Motor Deficits in Parkinson's Disease

**Code and resources for Fushiki et al. (2024), _"A Vulnerable Subtype of Dopaminergic Neurons Drives Early Motor Deficits in Parkinson's Disease."_**

[![Preprint](https://img.shields.io/badge/bioRxiv-v2%20%7C%2010.1101%2F2024.12.20.629776-b31b1b.svg)](https://www.biorxiv.org/content/10.1101/2024.12.20.629776v2)
[![Data: GEO](https://img.shields.io/badge/GEO-GSE285508-blue.svg)](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285508)
[![Data: Zenodo](https://img.shields.io/badge/Zenodo-20766032-1682D4.svg)](https://zenodo.org/records/20766032)
[![Mice: JAX](https://img.shields.io/badge/JAX-040185%20%7C%20040184-006747.svg)](https://www.jax.org/strain/040185)
[![License: CC BY 4.0](https://img.shields.io/badge/License-CC%20BY%204.0-lightgrey.svg)](http://creativecommons.org/licenses/by/4.0/)

---

## Overview

This repository contains the analysis code and pointers to all associated data and experimental
protocols for our study of selective dopaminergic neuron vulnerability in Parkinson's disease (PD).

Using MitoPark mice as a progressive PD model, we identify an *Anxa1*+ subtype of *Sox6*+
dopaminergic neurons (DANs) in the ventral substantia nigra pars compacta that is lost earlier than
other DAN populations, on a timescale matching the onset of motor symptoms. We generated an
*Anxa1*-Cre knock-in line to map the inputs and outputs of this population, and showed that
suppressing transmitter release from *Anxa1*+ neurons is sufficient to produce bradykinesia and
tremor — hallmarks of early PD.

The code here covers three analysis streams: **behavior**, **connectivity/anatomy**, and
**single-nucleus RNA sequencing**.

📄 **Preprint (current version — v2, posted 17 February 2026):**
<https://www.biorxiv.org/content/10.1101/2024.12.20.629776v2>
The DOI [10.1101/2024.12.20.629776](https://doi.org/10.1101/2024.12.20.629776) always resolves to the
latest version; [v1](https://www.biorxiv.org/content/10.1101/2024.12.20.629776v1) (20 December 2024)
remains available. The code in this repository corresponds to **v2**.

---

## Repository structure

```
FushikiA_etal_2024/
├── Behavioral_Analysis/     # Open field, cylinder, and Motion sensor analysis
├── Connectivity_Analysis/   # Anterograde/retrograde tracing, whole-brain registration & quantification
├── Sequencing_Analysis/     # snRNA-seq preprocessing, clustering, differential abundance (Milo)
├── LICENSE.txt
└── README.md
```

| Directory | Contents | Related figures |
| --- | --- | --- |
| [`Behavioral_Analysis/`](./Behavioral_Analysis) | Scripts for quantifying locomotion, rearing, tremor, and Pavlovian conditioning | _Figs. 1, 6_ |
| [`Connectivity_Analysis/`](./Connectivity_Analysis) | Registration to the Allen CCF and quantification of axonal projections and rabies-labeled inputs | _Fig. 2, 5_ |
| [`Sequencing_Analysis/`](./Sequencing_Analysis) | snRNA-seq QC, integration, subtype annotation, and differential abundance testing | _Figs. 2–4_ |

Figure numbers refer to the bioRxiv preprint (v2). Each analysis directory has its own `README.md`
listing the software environment, exact package versions, and the order in which to run the scripts.
Behavioral and connectivity analyses run in Python 3.11.5; sequencing analyses run in R 4.3.3.

---

## Data availability

| Resource | Repository | Accession |
| --- | --- | --- |
| Whole-brain imaging datasets | Brain Image Library | [10.35077/g.1209](https://doi.brainimagelibrary.org/doi/10.35077/g.1209) |
| Raw behavioral data and video recordings | BioStudies (EMBL-EBI) | [S-BSST2649](https://www.ebi.ac.uk/biostudies/studies/S-BSST2649) |
| Single-nucleus RNA-sequencing data | Gene Expression Omnibus | [GSE285508](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE285508) |
| Processed and supporting datasets | Zenodo | [20766032](https://zenodo.org/records/20766032) |

---

## Mouse lines

Both knock-in lines generated in the Costa lab and used in this study are openly available from
The Jackson Laboratory:

| Line | Full strain name | Allele | JAX stock | Description |
| --- | --- | --- | --- | --- |
| ***Anxa1*-Cre** (Anxa1-P2A-Cre) | C57BL/6-*Anxa1<sup>em1(cre)Rcos</sup>*/J | *Anxa1<sup>em1(cre)Rcos</sup>* | [040185](https://www.jax.org/strain/040185) | Cre recombinase expressed from the endogenous *Anxa1* promoter, giving genetic access to the *Anxa1*+ subtype of SNc dopaminergic neurons |
| **DAT-FlpO** (DAT-Flp, Slc6a3-T2A-FlpO) | C57BL/6-*Slc6a3<sup>em1(flpo)Rcos</sup>*/J | *Slc6a3<sup>em1(flpo)Rcos</sup>* | [040184](https://www.jax.org/strain/040184) | Codon-optimized FlpO expressed from the endogenous *Slc6a3* (DAT) promoter, labeling midbrain dopaminergic neurons |

Combined with the DAT-FlpO line, subtype-specific Cre drivers enable intersectional (Cre ∩ Flp)
targeting of molecularly defined dopaminergic subtypes — including *Anxa1*+ neurons via the
*Anxa1*-Cre line above, and other subtypes such as *Calb1*+ and *Vglut2*+ (*Slc17a6*+) neurons via
existing Cre lines. This is the strategy used for the tracing and behavioral experiments in this
study.

**Generation and validation.** The *Anxa1*-Cre line is described in this study. The DAT-FlpO line was
generated and validated alongside four other knock-in lines in:

> Albarran E, Fushiki A, Nelson A, Ng D, Chaimowitz C, Nikoobakht L, Sippy T, Peterka DS, Costa RM.
> *Generation of knock-in Cre and FlpO mouse lines for precise targeting of striatal projection
> neurons and dopaminergic neurons.* eLife Reviewed Preprint (2025).
> doi: [10.7554/eLife.108458](https://doi.org/10.7554/eLife.108458) ·
> [elifesciences.org/reviewed-preprints/108458](https://elifesciences.org/reviewed-preprints/108458)

---

## Experimental protocols

All wet-lab protocols are openly available on protocols.io:

- [Isolation of nuclei from adult mouse brain tissue (V.2)](https://www.protocols.io/view/isolation-of-nuclei-from-adult-mouse-brain-tissue-14egn8b36g5d/v2)
- [Costa Lab Nanoject virus injections — CU](https://www.protocols.io/view/costa-lab-nanoject-virus-injections-cu-kxygxy6j4l8j/v1)
- [Histology — Fushiki et al. 2024](https://www.protocols.io/view/histology-fushiki-et-al-2024-5qpvo9o49v4o/v1)
- [Open field & cylinder test behavior assays](https://www.protocols.io/view/open-field-amp-cylinder-test-behavior-assays-eq2ly6ypegx9/v1)
- [Pavlovian protocol — Fushiki et al. 2024](https://www.protocols.io/view/pavlovian-protocol-fushiki-et-al-2024-yxmvmdk79v3p/v1)

---

## Citation

If you use this code or data, please cite:

> Fushiki A, Ng D, Lewis ZR, Yadav A, Saraiva T, Hammond LA, Wirblich C, Tasic B, Menon V,
> Alves da Silva J, Costa RM. *A Vulnerable Subtype of Dopaminergic Neurons Drives Early Motor
> Deficits in Parkinson's Disease.* bioRxiv 2024.12.20.629776 [Preprint]. Version 2, posted
> 17 February 2026. doi: [10.1101/2024.12.20.629776](https://doi.org/10.1101/2024.12.20.629776)

<details>
<summary>BibTeX</summary>

```bibtex
@article{Fushiki2024,
  title   = {A Vulnerable Subtype of Dopaminergic Neurons Drives Early Motor Deficits in Parkinson's Disease},
  author  = {Fushiki, Akira and Ng, David and Lewis, Zachary R. and Yadav, Archana and
             Saraiva, Tatiana and Hammond, Luke A. and Wirblich, Christoph and Tasic, Bosiljka and
             Menon, Vilas and Alves da Silva, Joaquim and Costa, Rui M.},
  journal = {bioRxiv},
  year    = {2024},
  doi     = {10.1101/2024.12.20.629776},
  note    = {Preprint, version 2, posted 17 February 2026},
  url     = {https://www.biorxiv.org/content/10.1101/2024.12.20.629776v2}
}
```

</details>

---

## Contact

For questions about the code or data, please
[open an issue](https://github.com/afushiki/FushikiA_etal_2024/issues).

---

## Acknowledgements

This work was supported by grants from Aligning Science Across Parkinson's (ASAP-020551) through
the Michael J. Fox Foundation for Parkinson's Research (MJFF), the Parkinson's Disease Foundation
(PDFPF-RCE-1948), the National Institutes of Health (NIH, U19: 5U19NS104649-03), and a NARSAD Young
Investigator Grant from the Brain & Behavior Research Foundation (30086). This research was funded
in part through the NIH/NCI Cancer Center Support Grant P30CA013696 and used the Genomics and High
Throughput Screening Shared Resource. Additional support was received from the National Center for
Advancing Translational Sciences, NIH, under Grant Number UL1TR001873. Rabies viruses were produced
by the Center for Neuroanatomy with Neurotropic Viruses (CNNV), supported by the NIH under Grant
P40 OD010996.

For the purpose of open access, the authors have applied a CC BY public copyright license to all
Author Accepted Manuscript versions arising from this submission.

---

## License

[![CC BY 4.0](https://i.creativecommons.org/l/by/4.0/88x31.png)](http://creativecommons.org/licenses/by/4.0/)

This work is licensed under a
[Creative Commons Attribution 4.0 International License](http://creativecommons.org/licenses/by/4.0/).

---

## Changelog

| Date | Change |
| --- | --- |
| August 2026 | Added data, protocol, and mouse line links; revised README (A.F.) |
| February 2026 | Preprint updated to version 2 on bioRxiv (A.F.) |
| January 2026 | Added analysis scripts for the manuscript (A.F.) |
| January 2025 | Repository created (A.F.) |
