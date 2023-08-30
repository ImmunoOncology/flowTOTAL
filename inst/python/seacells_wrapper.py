# seacells_wrapper.py
import scanpy as sc
import numpy as np
import pandas as pd
import time
import matplotlib.pyplot as plt
import seaborn as sns
import anndata as ann
from scipy.sparse import csr_matrix
from collections import defaultdict
import SEACells
from SEACells.core import SEACells


def run_min_max_sampling_wrapper(file_counts, output, n_SEACells, n_waypoint_eigs=5, n_comps_pca=10, build_kernel_on="X_umap"):
  counts = csr_matrix(csr_matrix(np.loadtxt(file_counts), dtype=np.float32), dtype=np.float32)
  ad = ann.AnnData(counts)

  sc.pp.normalize_per_cell(ad)
  sc.pp.log1p(ad)
  sc.pp.neighbors(ad)
  sc.tl.umap(ad)
  sc.tl.pca(ad, n_comps=n_comps_pca, use_highly_variable=False)

  from SEACells.core import SEACells
  model = SEACells(ad, build_kernel_on=build_kernel_on, n_SEACells=n_SEACells, n_waypoint_eigs=n_waypoint_eigs,convergence_epsilon = 1e-5)
  model.construct_kernel_matrix()
  M = model.kernel_matrix
  model.initialize_archetypes()

  file_name = output + "/min-max-sampling-index.txt"

  with open(file_name, 'w') as file:
    file.write('\n'.join(str(cell) for cell in model.archetypes))


  return model.archetypes