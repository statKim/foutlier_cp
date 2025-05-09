# Conformal Outlier Detection for Multivariate Functional Data

This is the R code to implement the simulations and real data analysis of the following paper:
> Hyunsung Kim and Junyong Park (2025+). Conformal Outlier Detection for Multivariate Functional Data, *manuscript*.

## Description

- **R/foutlier_cp.R**: R functions to implement the proposed method and other functions are included.
- **py/adhd200_preprocessing.ipynb**: A notebook file for preprocessing the ADHD-200 dataset using `Nilearn` package in Python.
- **adhd_200_preprocess.R**: R code for converting the data structure of preprocessed ADHD-200 dataset from `Nileran`.
- **sim_clean.R**: R code of the simulation for clean training set (test outlier detection)
- **sim_mixed.R**: R code of the simulation for mixed training set (test outlier detection)
- **sim_mixed_train.R**: R code of the simulation for mixed training set (training outlier detection)
- **real_adhd200.R**: R code of the real data analysis


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
