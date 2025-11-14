# setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)

# Load data
seurat_integrated <- readRDS("14_INF_integration/out/INF_seurat_integrated_v5_RNA.rds")

# Settings
DefaultAssay(seurat_integrated) <- "RNA"

Reductions(seurat_integrated)

reduction <- "RNA_umap.harmony"
cluster.name <- "RNA_harmony_clusters_res.1"
# seurat_integrated$RNA_harmony_clusters_res.1

# Final UMAP
DimPlot(seurat_integrated, reduction = reduction, group.by = cluster.name, label = TRUE) +
  labs(title = glue("UMAP - post integration"),
       subtitle = cluster.name)

ggsave(glue(glue("16_INF_topDEGs/plot/UMAP_{cluster.name}.pdf")), 
       width = 8, 
       height = 7)

# Pre annotation 
DimPlot(seurat_integrated, reduction = reduction, group.by = "annotation", label = TRUE) +
  labs(title = glue("UMAP - post integration"),
       subtitle = "Annotation from study")

ggsave(glue(glue("16_INF_topDEGs/plot/UMAP_sample.pdf")), 
       width = 8, 
       height = 7)

# Join layers 
seurat_integrated[["RNA"]] <- JoinLayers(seurat_integrated[["RNA"]])
Layers(seurat_integrated[["RNA"]])

# Find markers 
# markers_5_vs_all <- FindMarkers(seurat_integrated, 
#                                 ident.1 = 5,
#                                 group.by = cluster.name)
# head(markers_5_vs_all, n = 100)

# Find all markers 
all_markers <- FindAllMarkers(seurat_integrated, 
                              group.by = cluster.name)

# # Ensure the gene names are in a column (if they were in row names)
df_all_markers <- all_markers #%>%
  # rownames_to_column("gene")

# Process markers
top_markers_list <- df_all_markers %>%
  filter(p_val_adj < 0.05) %>% # Filter for significant markers
  group_by(cluster) %>%
  # Sort by adjusted p-value (most significant) and then avg_log2FC (highest expression)
  arrange(p_val_adj, desc(avg_log2FC), .by_group = TRUE) %>%
  slice_head(n = 100) %>%
  ungroup() %>%
  # The split() function creates the named list needed for openxlsx
  split(., .$cluster)

# Export xlsx file 
out_file <- "16_INF_topDEGs/out/INF_Top_100_Cluster_DE_Markers.xlsx"

# Use openxlsx::write.xlsx, which takes the named list and writes
# each element as a separate sheet (sheet name = list name, i.e., Cluster ID)
openxlsx::write.xlsx(
  x = top_markers_list,
  file = out_file,
  overwrite = TRUE # Overwrite the file if it already exists
)



