# CNAttack
Repository for CSCB Final Project, Group 6

## Description
CNAttack is a computational tool to identify CNVs. This package is a CNA detection choose your own adventure that allows users to cutomize the package as follows: 

Choose your own clustering variable
* Cell Type
* Time point
* Conditioned group
  
Choose your own expression average
* Cluster specific
* Reference specific
* Global


## Installation
1. Clone the respository:
  git clone https://github.com/pfkrajca/CNAttack.git
3. Navigate to the project directory:
   cd CNAttack
5. Install dependencies:
   pip install TBD

## Usage
1. Import the package:
   from CNAttack import TBD
2. Use the package:
   result = TBD.TBD(arg1, arg2)
   print(result)

## Setup

Load the data.
```
adata = sc.read_h5ad("PBMC_simulated_cnas_041025.h5ad")
```

Before running the package ensure your adata object contains:

* `adata.var` must contain the columns `"gene_name"`, `"chromosome"`, `"start"`, and `"end"`.

* `adata.obs` must contain a clustering variable column such as `"cell_type"` that classifies each cell.


Import the required packages.
```
import os
import sys
import argparse
import numpy as np
import pandas as pd
import scanpy as sc
import scFates as scf
from anndata import AnnData
from typing import Dict, Tuple
from hmmlearn import hmm
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
import warnings
import matplotlib.pyplot as plt
from scipy.sparse import issparse
from scipy.stats import zscore
```


## Finding a Reference Cluster
This function plots the average gene expression across genomic positions for each cell type. These plots will help you visually identify which cell type has the lowest overall expression across the genome, which will serve as the reference cluster in later steps. You can input the name of your clustering variable into `group_by`. Default for all `group_by` is `"cell_type"`

```
plot_cluster_means_by_genomic_position(adata, layer="counts", group_by="cell_type", downsample=100)

```
Based off these plots determine a cell type to set as the baseline. Reference cluster should have lowest overall gene expression across the entire genome.


## Gene Expression Smoothing
After identifying the reference cluster, compute the smoothed expression profiles. Replace "reference_cluster" with the name of the cluster identified in the previous step. It returns the smoothed expression data, the average expression by cell type, the global average expression across all cells, and the average for the reference cluster.
```
smoothed, celltype_avg, global_avg, cluster_avg = compute_smoothed_profiles_from_adata(adata, group_by="cell_type", ref="reference_cluster")
```

Compute z-scores to measure deviations in expression for each cell relative to the smoothed reference profiles. It stores the resulting z-score matrices in adata.obsm, such as `adata.obsm["delta_global_z"]`.
```
adata= compute_all_cell_zscores_to_adata_optimized(smoothed, celltype_avg, global_avg, cluster_avg, adata, group_by="cell_type")
```

Filter z-scores and count significant values. Will need to run function three times for each reference profile z score matrices.
```
filtered, up, down = filter_and_count_zscores(
    adata.obsm["delta_global_z"],
    upper_thresh=1.5,
    lower_thresh=-1.5,
    min_cells=600,
    adata=adata,
    obsm_key="filtered_global_z"
)

filtered, up, down = filter_and_count_zscores(
    adata.obsm["delta_celltype_z"],
    upper_thresh=1.5,
    lower_thresh=-1.5,
    min_cells=600,
    adata=adata,
    obsm_key="filtered_celltype_z"
)

filtered, up, down = filter_and_count_zscores(
    adata.obsm["delta_cluster_z"],
    upper_thresh=1.5,
    lower_thresh=-1.5,
    min_cells=600,
    adata=adata,
    obsm_key="filtered_cluster_z"
)
```


## HMM-Based CNV Detection

Infer CNVs using Hidden Markov Model. This function trains an HMM on the filtered z-score matrix ("filtered_cluster_z") and classifies genomic regions into one of three states: 0 for loss, 1 for neutral, and 2 for gain. You can adjust n_iter for more iterations, or chunk_size for memory-efficient processing. The inferred CNV states are stored in adata.var with column names prefixed by output_prefix ("hmm_cnv_state").  

```
adata = detect_cnvs_with_hmm_final(adata, matrix_name="filtered_cluster_z", n_components=3,
                             n_iter=50, random_state=42, output_prefix="hmm_cnv",
                             min_non_nan=20, chunk_size=100
)
```

Format the results. Function merges neighboring CNVs based on genomic proximity and filters them by size and z-score deviation. It returns a summary of CNV regions including genomic coordinates, affected genes, and cell counts. The CNV region metadata is saved in `adata.obs["detected_cnvs"]`

```
adata, cnv_stats = format_detected_cnvs_with_cell_counts(
    adata,
    state_col='hmm_cnv_state',
    output_col='detected_cnvs',
    group_by_col='cell_type'
    max_gap=500000,
    min_region_size=5000,
    z_threshold=1.5
)


```



