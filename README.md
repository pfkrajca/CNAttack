# group6-cscb
Repository for CSCB Final Project, Group 6

## Description

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

Before running the package ensure your adata object contains:

  `adata.var` with the following columns:

  *  `"gene_name"`

  *  `"chromosome"`

  *  `"start"`

  *  `"end"`

  `adata.obs` with a `"cell_type"` column.




## Finding a Reference Cluster

Import the required packages.
`import os, sys
os.environ['R_HOME'] = sys.exec_prefix+"/lib/R/"
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
from scipy.sparse import issparse`

## Gene Expression Smoothing

## HMM-Based CNV Detection


## Examples

## License
This project is licensed under the MIT License, see the LICENSE file for more details.

## Support
If you have any questions or concerns, you can open an issue or contact us at TBD


