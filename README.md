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

General CNAttack Pipeline
1. **Plotting**: Visualize gene expression across genomic positions
2. **Smoothing**: Apply smoothing windows to expression profiles
3. **Z-score calculation**: Compute expression deviations from expected values
4. **Filtering**: Focus on significant expression changes
5. **HMM detection**: Use Hidden Markov Models to identify CNV patterns
6. **Annotation**: Format, merge, and annotate detected CNVs

## Example Data Sets
Sample scRNA-seq data sets used in this repository can be found at this public Google Drive link: https://drive.google.com/drive/folders/198Vcme8ewUsUOdKq9GVszRskIuurNxQT?usp=sharing 

## Installation

A usage example of CNAttack for command-line interfaces is provided in example_CNA.py. Please note that we generated and tested our code in Google Colab, therefore, there may be some differences in the conversion of our package from Google Colab to Argparse formatting for use in the command line. The below instructions are more attuned to Google Colab users. Instructions for the example_CNA.py file, can be found after this.

1. **Clone the repository**:
    ```bash
    git clone https://github.com/pfkrajca/CNAttack.git
    ```

2. **Navigate to the project directory**:
    ```bash
    cd CNAttack
    ```

3. **Install dependencies**:
    ```bash
    pip install -r requirements.txt
    ```

4. **Import the package in your Python code**:
    ```python
    from CNAttack import filter_zscores
    from CNAttack import format_CNVs
    from CNAttack import hmm_for_CNAs
    from CNAttack import plot_cluster
    from CNAttack import smooth_profiles
    from CNAttack import zscores_compute
    from CNAttack.plot_cluster import plot_cluster_means_by_genomic_position
    from CNAttack.smooth_profiles import compute_smoothed_profiles_from_adata
    from CNAttack.zscores_compute import compute_all_cell_zscores_to_adata_optimized
    from CNAttack.filter_zscores import filter_and_count_zscores
    from CNAttack.hmm_for_CNAs import detect_cnvs_with_hmm_final
    from CNAttack.format_CNVs import format_detected_cnvs_with_cell_counts
    ```
    
## Setup & Usage

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
## Command Line Interface Example (example_CNAttack.py)

## Quick Start

The simplest way to run the complete workflow is using the example script:

```bash
python example_CNAttack.py input.h5ad --output results/my_analysis
```
The unput .h5ad file should be the PBMC simulated data set provided. This can be found ___. 
This way, the same annData file will be use throughout the whole process to ensure clean detection of CNAs.

## Required Data Format

Your input file should be an AnnData object (.h5ad) with:
- Gene expression data in the main matrix or specified layer, like 'counts'
- Gene annotations in `adata.var` including:
  - `chromosome`: Chromosome location (e.g., "1", "X", "chrX")
  - `start`: Gene start position (integer)
  - `end`: Gene end position (integer)
- Cell annotations in `adata.obs` including:
  - Cell type or grouping information (default: "cell_type" column)
  - Other possible clustering parameters could include "time_point" or "condition"

## Example Workflow

```python
# Complete workflow
python cnattack_example.py input.h5ad --output results/full_analysis --plot

# Run with custom parameters
python cnattack_example.py input.h5ad \
    --output results/custom_analysis \
    --layer normalized \
    --group_by leiden_clusters \
    --window_size 50 \
    --ref_group normal_cells \
    --upper_thresh 1.5 \
    --lower_thresh -1.5 \
    --min_cells 300 \
    --n_components 3 \
    --z_threshold 1.8

# Run specific steps only
python cnattack_example.py input.h5ad \
    --output results/partial \
    --start_step 4 \
    --end_step 6
```

## Individual Modules

Each step can also be run individually using the module-specific scripts:

```bash
#Plot expression by genomic position
python -m CNAttack.plot_cluster input.h5ad --layer counts --group_by cell_type --downsample 100

#Format CNVs from an h5ad with HMM states
python -m CNAttack.format_CNVs hmm_results.h5ad --output_col detected_cnvs
```

## Output Files

The workflow generates several files with the specified prefix:

- `*_smoothed.h5ad`: AnnData with smoothed expression profiles
- `*_zscores.h5ad`: AnnData with computed z-scores
- `*_filtered.h5ad`: AnnData with filtered significant genes
- `*_hmm.h5ad`: AnnData with HMM state assignments
- `*_cnvs.h5ad`: Final AnnData with CNV annotations
- `*_cnv_regions.csv`: Table of detected CNV regions with statistics

In Google Colab, it is easier to understand the processing of your data than in the command line, the integration of file generation for the command line specific example ensures the user can check their work. This is critically important as our package is a choose your own adventure!
