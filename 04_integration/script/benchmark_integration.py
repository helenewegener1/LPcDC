#!/usr/bin/env python3

# https://scib-metrics.readthedocs.io/en/stable/notebooks/lung_example.html

import scanpy as sc
import pandas as pd
from scib_metrics.benchmark import Benchmarker

# 1. Define input data
files = {
    "RNA_cca": "/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/mydata_RNA_integrated.cca.h5ad",
    "RNA_mnn": "/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/mydata_RNA_integrated.mnn.h5ad",
    "RNA_harmony": "/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/mydata_RNA_integrated.harmony.h5ad",
    "SCT_mnn": "/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/mydata_SCT_integrated.mnn.h5ad",
    "SCT_harmony": "/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/mydata_SCT_integrated.harmony.h5ad",
}

# 2. Load all datasets into a dictionary
adatas = {}
for name, path in files.items():
    adata = sc.read_h5ad(path)
    # Find the integrated embedding and rename it to a unique key
    for key in list(adata.obsm.keys()):
        print(key)
        adata.obsm[f"X_{name}"] = adata.obsm.pop(key)
    adatas[name] = adata

print('--------------step 4---------------')

for name, adata in adatas.items():
    print(name, adata.obsm.keys())

print('----------------------------------')

print(adatas["SCT_mnn"].obsm.keys())

# # 3. Define batch and label columns
# # Make sure these columns exist in your Seurat object's metadata before export!
# batch_key = "orig.ident"   # adjust if your batch column has a different name
# # label_key = "celltype"     # use ground-truth celltype annotation - if any

# # Create a dummy variable 
# for adata in adatas.values():
#     adata.obs["celltype"] = ["Unknown"] * adata.n_obs  # Create a dummy label column

# # 4. Initialize Benchmarker
# bm = Benchmarker(
#     adata,
#     batch_key = batch_key,
#     label_key="cell_type",
#     embedding_obsm_keys=["X_RNA_cca", "X_RNA_mnn", "X_RNA_harmony", "X_SCT_mnn", "X_SCT_harmony"],
#     n_jobs=6,
# )
# bm.benchmark()

# # Now pass embeddings when calling benchmark()
# bm.benchmark(embeddings=[f"X_{name}" for name in adatas.keys()])

# # Collect and save results
# df = bm.get_results(min_max_scale=False)
# print(df)
# df.to_csv("/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/out/integration_benchmark_results.csv")

# # Save the results plot
# fig = bm.plot_results_table()
# fig.savefig("/Users/srz223/Documents/projects/project_cDC/LPcDC/04_integration/plot/integration_benchmark_results.png", dpi=300, bbox_inches="tight")

