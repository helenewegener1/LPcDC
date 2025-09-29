#!/usr/bin/env python3

# https://scib-metrics.readthedocs.io/en/stable/notebooks/lung_example.html

import numpy as np
import scanpy as sc
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection


# 1. Define input data
path = "/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/mydata_RNA.h5ad"

# 2. Load all datasets into a dictionary
adata = sc.read(path)
print(adata)

print('----------------------------------')

# Renaming
adata.obsm["Unintegrated"] = adata.obsm["X_pca"]
print(adata)

print('--------------step 4---------------')

# Perform the benchmark
bm = Benchmarker(
    adata,
    batch_key="orig.ident",
    label_key="cell_type",
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    embedding_obsm_keys=["X_integrated.cca", "X_integrated.harmony", "X_integrated.mnn", "Unintegrated"],
    n_jobs=6,
)
bm.benchmark()


print('----------------------------------')



# Collect and save results
df = bm.get_results(min_max_scale=True)
print(df)
df.to_csv("/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/integration_benchmark_results_min_max_scale.csv")

df = bm.get_results(min_max_scale=False)
print(df)
df.to_csv("/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/integration_benchmark_results.csv")

# Save the results plot
# import matplotlib.pyplot as plt

bm.plot_results_table(min_max_scale=True, save_dir="/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/plot/integration_benchmark_results_min_max_scale.png")

bm.plot_results_table(min_max_scale=False, save_dir="/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/plot/integration_benchmark_results.png")
