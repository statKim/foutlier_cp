# Conformal Outlier Detection for Multivariate Functional Data

This is the R code to implement the simulations and real data analysis of the following paper:
> Hyunsung Kim and Junyong Park (2025+). Conformal Outlier Detection for Multivarate Functional Data, *submitted*.

## Description

- **R/foutlier_cp.R**: R functions to implement the proposed method and other functions are included.
- **R/other_functions.R**: functions for simulation data, performance measure, and etc.
- **py/adhd200_preprocessing.ipynb**: A notebook file for preprocessing the ADHD-200 dataset using `Nilearn` package in Python.

Additionally, you need to download the modified version of the R package `mrfDepth`:

```R
devtools::install_github("statKim/mrfDepth")
```

## Data source

- **ADHD-200**
  - Download from NITRC website: https://www.nitrc.org/frs/?group_id=383
  - Descriptions of preprocessing: https://www.nitrc.org/plugins/mwiki/index.php?title=neurobureau:AthenaPipeline
  - We downloaded the following preprocessed data.
    - Peking_1_preproc_filtfix.tar
    - Peking_2_preproc_filtfix.tar
    - Peking_3_preproc_filtfix.tar
