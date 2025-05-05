# group6-cscb
Repository for CSCB Final Project, Group 6

## Description
___ is a computational tool to identify CNVs. 

## Installation
1. Clone the respository:
  git clone https://github.com/pfkrajca/group6-cscb.git
3. Navigate to the project directory:
   cd group6-cscb
5. Install dependencies:
   pip install TBD

## Usage
1. Import the package:
   from group6-cscb import TBD
2. Use the package:
   result = TBD.TBD(arg1, arg2)
   print(result)

## Setup

Load data.
```
adata = sc.read_h5ad("PBMC_simulated_cnas_041025.h5ad")
```

Before running the package ensure your adata object contains:

* `adata.var` must contain the columns `"gene_name"`, `"chromosome"`, `"start"`, and `"end"`.

* `adata.obs` must contain a `"cell_type"` column that identifies the cell type for each cell.


Import the required packages.
```
import os, sys
import scanpy as sc
import scFates as scf
import numpy as np
import pandas as pd
from anndata import AnnData
from typing import Dict, Tuple
from hmmlearn import hmm
from sklearn.preprocessing import StandardScaler
from tqdm import tqdm
import warnings
import matplotlib.pyplot as plt
from scipy.sparse import issparse
```


## Finding a Reference Cluster
This function plots the average gene expression across genomic positions for each cell type. The layer parameter specifies which data layer (e.g., "counts") to use. The downsample parameter reduces the number of cells used for plotting (default is 100) to improve performance for large datasets. These plots help you visually identify which cell type has the lowest overall expression across the genome, which will serve as the reference cluster in later steps.

```
plot_cluster_means_by_genomic_position(adata, layer="counts", downsample=100)

```
Based off these plots determine a cell type to set as the baseline. Reference cluster should have lowest overall gene expression across the entire genome.


## Gene Expression Smoothing
After identifying the reference cluster, compute the smoothed expression profiles. Replace "reference_cluster" with the name of the cluster identified in the previous step. It returns the smoothed expression data, the average expression by cell type, the global average expression across all cells, and the average for the reference cluster.
```
smoothed, celltype_avg, global_avg, cluster_avg = compute_smoothed_profiles_from_adata(adata, ref="reference_cluster")
```

Next, compute z-scores to measure deviations in expression for each cell relative to the smoothed reference profiles. This function calculates z-scores comparing each cell’s expression to the cell-type-specific average, the global average, and the reference cluster average. It stores the resulting z-score matrices in adata.obsm, such as `adata.obsm["delta_global_z"]`.
```
adata= compute_all_cell_zscores_to_adata_optimized(smoothed, celltype_avg, global_avg, cluster_avg, adata)
```





## HMM-Based CNV Detection





## Examples

## License
This project is licensed under the MIT License, see the LICENSE file for more details.

## Support
If you have any questions or concerns, you can open an issue or contact us at TBD


