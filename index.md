---
layout: page
title: Reproducibility Guide
---

## Table of Contents

- [I. Project Information](#i-project-information)
- [II. CUBIC Project Directory Structure](#ii-cubic-project-directory-structure)
- [III. Code Documentation](#iii-code-documentation)
  - [0. Install required packages and libraries](#0-install-required-packages-and-libraries)
  - [1. Unzip datalad cloned data directories and copy required files](#1-unzip-datalad-cloned-data-directories-and-copy-required-files)
  - [2. Calculate intrinsic timescale for fMRI time series](#2-calculate-intrinsic-timescale-for-fmri-time-series)
  - [3. Concatenate the timescale and QC data](#3-concatenate-the-timescale-and-qc-data)
  - [4. Sample selection](#4-sample-selection)
  - [5. Run age analysis using Generalized Additive Models (GAMs)](#5-run-age-analysis-using-generalized-additive-models-gams)
  - [6. Run additional analysis and visualize the results](#6-run-additional-analysis-and-visualize-the-results)

## I. Project Information

**Abstract**

Intrinsic timescale is a commonly used measure of spontaneous neural dynamics that quantifies the timescale of processing of neuronal populations. Intrinsic timscale displays a hierarchical organization across multiple species and imaging modalities, with shorter timescales in sensorimotor cortex compared to association cortex. However, less is known about how intrinsic timescale evolves during development in the human brain and whether its cortical maturation patterns generalize to independent developmental samples. Here we estimate the intrinsic timescale in multiple datasets of youth (HCPD: n=565; HBN: n=729; age range 8–22 years) and investigate its neurodevelopmental patterns. We find that developmental changes in the intrinsic timescale follow a hierarchical pattern that recapitulates an archetypical axis spanning sensorimotor to association cortices (S–A axis). Our analysis of an independent healthy young adult dataset (HCPYA: n=973, age range 22–37 years) underscores the specificity of the developmental findings, suggesting that the intrinsic timescale develops along the S–A axis in youth and stabilizes in adulthood. Together, these results reveal convergence between major axes of cortical organization and development, highlighting intrinsic timescale as a principled marker of hierarchical brain maturation in youth.


 
**Team**

| **Role** | **Name** |
| --- | --- |
| **Project Lead** | Golia Shafiei |
| **Faculty Lead** | Theodore D. Satterthwaite |
| **Analytic Replicator** | Joëlle Bagautdinova |
| **Collaborators** | Valerie J. Sydnor, Dani S. Bassett5, Deanna M. Barch, Matthew Cieslak, Yong Fan, Elizabeth Flook, Alexandre R Franco, Gregory Kiar, Audrey C. Luo, Michael P Milham, Linden Parkes, Taylor Salo, Leah H. Somerville, Tien T. Tong, Russell T. Shinohara |

**Project Timeline**

| **Project Start Date** | November 2023 |
| --- | --- |
| **Current Project Status** | [Preprint on bioRxiv](https://doi.org/10.64898/2026.04.08.717312) |

**Code and Communication**

| **Github Repository** | [https://github.com/PennLINC/shafiei_timescale/tree/main](https://github.com/PennLINC/shafiei_timescale/tree/main) |
| --- | --- |
| **Slack Channel** | #shafiei_timescales |

**Datasets**

- HCP-D as discovery developmental sample (replication); HBN as replication discover sample; HCP-YA as an independent sample for sensitivity analysis

<!-- **Conference Presentations**

-  -->


---

## II. CUBIC Project Directory Structure

The project directory on CUBIC is `/cbica/projects/developmental_gradients`. 

- Code for the final manuscript is in `/cbica/projects/developmental_gradients/gitrepo/shafiei_timescale/code/`
- The directory on Github is `~/shafiei_timescale/code/`

| **Directory** | **Description** |
| --- | --- |
| `~/shafiei_timescale/code/` | Main manuscript code |
| `~/shafiei_timescale/data/` | Note that this directory is not shared publicly given that it includes restricted behavioral and demographic data. It includes data derivatives for each dataset. Specifically, it includes estimated intrinsic timescale for each individual's fMRI time-series data. It also includes related demographics and quality control (QC) data. |
| `~/shafiei_timescale/results/` | This directory includes all the results (i.e., figures and CSV files). Note that this directory is not shared publicly either given the large file sizes. But all results and figures can be regenerated using the code. |

---

## III. Code Documentation
***Overview of Analytic Workflow***

| **Step** | **Task** | **Notes** |
| --- | --- | --- |
| 0 | Install required packages and libraries |  |
| 1 | Unzip datalad cloned data directories and copy required files (Python) |  |
| 2 | Calculate intrinsic timescale for fMRI time series (Python) |  |
| 3 | Concatenate the timescale and QC data (Python) |  |
| 4 | Sample selection (Python) |  |
| 5 | Run age analysis using Generalized Additive Models (GAMs) (R) |  |
| 6 | Run additional analysis and visualize the results (Python) |  |

### 0. Install required packages and libraries

- `~/shafiei_timescale/code/python_env` folder contains the required `.yml` file and a series of commands to create and activate the Python environment used in this project.
- `~/shafiei_timescale/code/r_env` folder contains the instructions and files required to create and activate the R environment used in this project. Note that the R environment installation was done locally on a MacOS.

### 1. Unzip datalad cloned data directories and copy required files

- After cloning and getting the preprocessed data using datalad, the Python script `scp_unzip_xcpdfiles_xcpd.py` is used to unzip the indiviudal-level data and copy the desired files for downstream analyses (i.e., `.ptseries.nii` and QC files).

### 2. Calculate intrinsic timescale for fMRI time series

- Once all the time-series data are copied over, `scp_timescale_acf.py` is used to calculate the intrinsic timescale using the autocorrelation function of fMRI time series for each indiviudal and brain region. Note this script takes some time to run. I would recommend running it via `screen`.

### 3. Concatenate the timescale and QC data

- `scp_concatenate_timescale.py` is then used to concatenate the estimated timescale data and the associated QC information across indiviudals for each dataset.

### 4. Sample selection

- `scp_sample_<dataset>.py` is then used to perform the sample selection procedure for each dataset. The sample selection scripts output `.tsv` files that are ready for age analysis in R.

### 5. Run age analysis using Generalized Additive Models (GAMs)

- The age analysis with GAMs is performed using R (version `4.2.2`). `scp_timescale.R` is the main script that runs GAMs for each dataset and generates a subset of final figures. `scp_timescale.R` calls the functions stored in `fcn_GAM_timescale.R`.

### 6. Run additional analysis and visualize the results

- A subset of final figures were generated during the age analysis in the previous step. `scp_timescale_plotting.py` is used to perform a series of additional analyses and generate the remaining figures.
