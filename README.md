# flowTOTAL: flow cyTometry auTOmaTic AnaLysis <a href='https://github.com/ImmunoOncology/flowTOTAL'><img src='man/figures/logo.png' align="right" height="139" /></a>

## Overview

flowTOTAL is a user-friendly command line workflow to analyze flow cytometry data. Preprocessing, conventional analysis, and unsupervised analysis are the three primary portions of the pipeline. The user must specify the folder containing the FCS files, the metadata associated with each file, and the back-gating marker as input. During preprocessing, each.FCS will be corrected for fluorescence spillover (compensation), abnormalities will be detected by evaluating flow rate and signal acquisition, and doublets will be removed using forward scatter (QC). A population of interest will be specified for the traditional analysis, and the pipeline will proceed with normalization, dimensionality reduction, and clustering.

## Installation

flowTOTAL will be available soon!!

Dummy:

```{R}

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

BiocManager::install("flowCore")
BiocManager::install("ggcyto")
BiocManager::install("openCyto")
BiocManager::install("flowAI")
BiocManager::install("PeacoQC")

remotes::install_github(repo = "ImmunoOncology/flowTOTAL", ref = "dev-1.0.0")

```

## Version 

0.02
