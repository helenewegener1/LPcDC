#!/usr/bin/env python3

# https://scib-metrics.readthedocs.io/en/stable/notebooks/lung_example.html

import numpy as np
import scanpy as sc
from scib_metrics.benchmark import Benchmarker, BioConservation, BatchCorrection


# 1. Define input data
# path = "/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/mydata_all.h5ad"
# path = "/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/mydata_v4.h5ad"
# path = "/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/mydata_RNA_v5.h5ad"
path = "/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/mydata_SCT_v5.h5ad"

# 2. Load all datasets into a dictionary
adata = sc.read(path)
print(adata)

print('----------------------------------')

print(adata.obsm)

# Renaming
adata.obsm["Unintegrated"] = adata.obsm["X_PCA"]
print(adata)

print('--------------step 4---------------')

# Perform the benchmark
bm = Benchmarker(
    adata,
    batch_key="orig.ident",
    label_key="seurat_clusters",
    bio_conservation_metrics=BioConservation(),
    batch_correction_metrics=BatchCorrection(),
    # embedding_obsm_keys=["Unintegrated", "X_RNA_integrated_cca", "X_RNA_integrated_harmony", "X_RNA_integrated_mnn", "X_RNA_integrated_rpca"],
    embedding_obsm_keys=["X_SCT_integrated_harmony", "X_SCT_integrated_cca", "X_SCT_integrated_rpca", "Unintegrated"],
    n_jobs=6,
)

bm.benchmark()


print('----------------------------------')



# Collect and save results
# df = bm.get_results(min_max_scale=True)
# print(df)
# df.to_csv("/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/integration_benchmark_results_min_max_scale.csv")

df = bm.get_results(min_max_scale=False)
print(df)
df.to_csv("/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/integration_benchmark_results.csv")

# Save the results plot
# import matplotlib.pyplot as plt

# bm.plot_results_table(min_max_scale=True, save_dir="/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/plot/integration_benchmark_results_min_max_scale.png")

bm.plot_results_table(min_max_scale=False, save_dir="/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/plot/integration_benchmark_results.png")
