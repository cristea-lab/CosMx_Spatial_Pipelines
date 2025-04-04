# CosMx_Spatial_Pipelines

**CosMx_Spatial_Pipelines** is a collection of R scripts designed to process and analyze spatial transcriptomics data, specifically using the CosMx 1000 genes panel from Nanostring. This repository includes pipelines for data quality control (QC), cell subtyping, and ongoing development for analyzing tumor polarity and fibroblast activity.

## Table of Contents

1. [Introduction](#introduction)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Scripts Overview](#scripts-overview)
    - [1. QC Preprocessing Single Sample Pipeline](#1-qc-preprocessing-single-sample-pipeline)
    - [2. InSituType Single Sample Processing Pipeline](#2-insitutype-single-sample-processing-pipeline)
    - [3. Tumor Polarity and Fibroblast Single Sample Processing Pipeline](#3-tumor-polarity-and-fibroblast-single-sample-processing-pipeline)
5. [Authors and Contributors](#authors-and-contributors)
6. [License](#license)

## Introduction

The **CosMx_Spatial_Pipelines** repository provides tools to analyze spatial single-cell RNA sequencing data. This set of scripts is designed to streamline data processing for cancer research using the CosMx 1000 genes panel, focusing on QC, subtyping, and in-depth spatial analysis.

## Installation

To run these scripts, you need to have R installed along with the necessary packages. You can install the required R packages using the following commands:

```R
install.packages(c("Seurat", "dplyr", "tidyverse", "here"))
BiocManager::install(c("sceasy", "reticulate", "SpatialDecon", "InSituType", "dittoSeq"))
```

**Required R Packages:**
- **Seurat**: For single-cell RNA-seq data analysis.
- **sceasy**: For converting between different single-cell data formats.
- **reticulate**: For interfacing with Python.
- **dplyr**: For data manipulation.
- **SpatialDecon**: For spatial deconvolution analysis.
- **InSituType**: For cancer cell subtyping.
- **tidyverse**: A collection of R packages designed for data science.
- **here**: For constructing paths to files.
- **dittoSeq**: For visualization and analysis of single-cell RNA-seq data.

## Usage

Each script is designed to be run from the command line with specific arguments. Below is an overview of each script and its usage.

## Scripts Overview

### 1. QC Preprocessing Single Sample Pipeline

**File:** `1_QC_preprocessing_single_sample_pipeline.R`

**Description:**  
This script performs quality control preprocessing on single-sample data obtained from AtoMx. It generates QC plots to help users determine appropriate QC thresholds and outputs a Seurat object of the sample.

**Usage:**  
```bash
Rscript 1_QC_preprocessing_single_sample_pipeline.R <sample_name> <input_folder> <output_folder>
```

#### Input Parameters

| Parameter       | Description                                    |
|-----------------|------------------------------------------------|
| `<sample_name>` | Name of the sample being processed.            |
| `<input_folder>`| Directory containing input data from AtoMx.    |
| `<output_folder>`| Directory where outputs will be stored.       |

#### Outputs

| File Name                              | Description                                                |
|----------------------------------------|------------------------------------------------------------|
| `1_1_QC_probe_count_histogram.png`     | Histogram of probe counts.                                 |
| `1_2_QC_feature_count_histogram.png`   | Histogram of feature counts.                               |
| `1_3_QC_negative_probe_count_histogram.png` | Histogram of negative probe counts.                    |
| `A_1_pre_QC_sample.rds`                | Seurat object of the sample pre-QC.                        |

### 2. InSituType Single Sample Processing Pipeline

**File:** `2_insitutype_single_sample_processing_pipeline.R`

**Description:**  
This script runs the InSituType algorithm from Nanostring for cancer cell subtyping based on the input sample data. It allows the user to specify the range of clusters for analysis and set QC thresholds.

**Usage:**  
```bash
Rscript 2_insitutype_single_sample_processing_pipeline.R <sample_name> <input_folder> <output_folder> <institute_profile_data> <cluster_range_low> <cluster_range_high> <min_nCount_RNA> <max_nFeature_negprobes>
```

#### Input Parameters

| Parameter                | Description                                                                                   |
|--------------------------|-----------------------------------------------------------------------------------------------|
| `<sample_name>`          | Name of the sample being processed.                                                           |
| `<input_folder>`         | Directory containing input data from AtoMx.                                                   |
| `<output_folder>`        | Directory where outputs will be stored.                                                       |
| `<institute_profile_data>` | File path to the InSituType profile external data from Nanostring.                      |
| `<cluster_range_low>`    | Lower end of the range for the number of clusters to run.                                      |
| `<cluster_range_high>`   | Higher end of the range for the number of clusters to run.                                     |
| `<min_nCount_RNA>`       | Minimum nCount RNA for QC threshold.                                                           |
| `<max_nFeature_negprobes>`| Maximum nFeature negative probes for QC threshold.                                          |

#### Outputs

- **Seurat Object:**  
  - `A_2_post_QC_and_normalization_20_counts_per_cell.rds`: Seurat object file after QC and normalization.

- **Cluster-Specific Folders:**  
  For each specified number of clusters, a folder (e.g., `/12_clusters`) is created. Inside each folder, files with corresponding results are generated:

  | File Name                                                        | Description                                                                                       |
  |------------------------------------------------------------------|---------------------------------------------------------------------------------------------------|
  | `2_1_flightpath_initial_cluster.png`                             | Initial cluster visualization plot.                                                               |
  | `2_2_flightpath_new_cluster.png`                                 | New cluster visualization plot after refinement.                                                  |
  | `2_3_dotplot_marker_genes.png`                                   | Dot plot of marker genes across clusters.                                                         |
  | `2_4_barplot_level_1_clusters.png`                               | Bar plot showing distribution of level 1 clusters.                                                |
  | `2_5_barplot_level_2_clusters.png`                               | Bar plot showing distribution of level 2 clusters.                                                |
  | `2_6_barplot_level_3_clusters.png`                               | Bar plot showing distribution of level 3 clusters.                                                |
  | `3_1_violin_plot_of_marker_scores_stratified_by_collapsed_cell_type_1.png` | Violin plot of marker scores stratified by collapsed cell type 1.                             |
  | `3_2_violin_plot_of_marker_scores_stratified_by_collapsed_cell_type_2.png` | Violin plot of marker scores stratified by collapsed cell type 2.                             |
  | `4_1_spatial_plot.png`                                           | Spatial plot of cell distribution.                                                                |
  | `4_2_spatial_plot_level_1_clusters.png`                          | Spatial plot showing level 1 cluster distribution.                                                |
  | `B_1_semi_sup_insitutype_refined.rds`                            | Refined InSituType results in an RDS file.                                                        |
  | `B_2_semi_sup_insitutype_fully_labeled.rds`                      | Fully labeled InSituType results in an RDS file.                                                  |

### 3. Tumor Polarity and Fibroblast Single Sample Processing Pipeline

**File:** `3_tumor_polarity_and_fibrolast_single_sample_processing_pipeline.R`

**Description:**  
This script is currently under development. It aims to analyze tumor polarity and fibroblast activity in single-sample spatial transcriptomics data.

**Usage:**  
Details will be provided upon completion of the script development.

## Authors and Contributors

- **Xi (Sam) Wang, SM**: Author and computational biologist at the Department of Data Science, Dana-Farber Cancer Institute, supervised by Dr. Simona Cristea, PhD.
- **Len Taing: Maintainer and scientific programmer from the Center for Functional Cancer Epigenetics, at Dana-Farber Cancer Institute.
- **Dr. Alexander Jordan, MD**: Contributor from the lab of Dr. Jonathan Nowak, MD, PhD, at Dana-Farber Cancer Institute.
- **Dr. Simona Cristea, PhD**: Contributor and supervisor, Department of Data Science, Dana-Farber Cancer Institute.
- **Sudheshna Bodapati**: Previous member whose work contributed to the development of these pipelines.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for more details.
