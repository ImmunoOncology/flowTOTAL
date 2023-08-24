[![Github All Releases](https://img.shields.io/github/downloads/ImmunoOncology/flowTOTAL/total.svg)]()


# flowTOTAL: flow cyTometry auTOmaTic AnaLysis <a href='https://github.com/ImmunoOncology/flowTOTAL'><img src='man/figures/logo.png' align="right" height="139" /></a>

## Overview

flowTOTAL is a user-friendly command line workflow to analyze flow cytometry data. Preprocessing, conventional analysis, and unsupervised analysis are the three primary portions of the pipeline. The user must specify the folder containing the FCS files, the metadata associated with each file, and the back-gating marker as input. During preprocessing, each.FCS will be corrected for fluorescence spillover (compensation), abnormalities will be detected by evaluating flow rate and signal acquisition, and doublets will be removed using forward scatter (QC). A population of interest will be specified for the traditional analysis, and the pipeline will proceed with normalization, dimensionality reduction, and clustering.

## Installation

flowTOTAL will be available soon!!

Dummy:

```{R}

if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")

devtools::install_github("RGLab/RProtoBufLib")
devtools::install_github("RGLab/cytolib")
devtools::install_github("RGLab/flowCore")
devtools::install_github("RGLab/flowWorkspace")
devtools::install_github("RGLab/ggcyto")
devtools::install_github("RGLab/openCyto")
devtools::install_github("RGLab/flowStats")
devtools::install_github("giannimonaco/flowAI")
devtools::install_github("saeyslab/PeacoQC")

remotes::install_github(repo = "ImmunoOncology/flowTOTAL", ref = "dev-1.0.0")

```

## Version 

0.02
