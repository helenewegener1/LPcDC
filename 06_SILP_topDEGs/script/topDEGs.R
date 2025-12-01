setwd("~/Documents/projects/project_cDC/LPcDC/")

# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)

# Load data
seurat_integrated <- readRDS("04_SILP_integration/out/SILP_seurat_integrated_v5_RNA.rds")

n_cell <- ncol(seurat_integrated)

# Settings
DefaultAssay(seurat_integrated) <- "RNA"

Reductions(seurat_integrated)

reduction <- "RNA_umap.harmony"
# cluster.name <- "RNA_harmony_clusters_res.0.8"
cluster.name <- "RNA_harmony_clusters_res.1"

# reduction <- "RNA_umap.cca"
# cluster.name <- "RNA_cca_clusters_res.0.8"
# seurat_integrated$RNA_cca_clusters_res.0.8

# Final UMAP
DimPlot(seurat_integrated, reduction = reduction, group.by = cluster.name, label = TRUE) +
  labs(title = glue("UMAP - post integration"),
       subtitle = cluster.name,
       caption = glue("N cell: {n_cell}")
       )

# FeaturePlot(seurat_integrated, features = "Prdm16", reduction = reduction)

ggsave(glue(glue("06_SILP_topDEGs/plot/UMAP_{cluster.name}.pdf")), 
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
out_file <- glue("06_SILP_topDEGs/out/SILP_Top_100_Cluster_DE_Markers_{cluster.name}.xlsx")

# Use openxlsx::write.xlsx, which takes the named list and writes
# each element as a separate sheet (sheet name = list name, i.e., Cluster ID)
openxlsx::write.xlsx(
  x = top_markers_list,
  file = out_file,
  overwrite = TRUE # Overwrite the file if it already exists
)

########################### Cell cycle regressed out ########################### 

seurat_object <- readRDS("05_SILP_annotation/out/seurat_integrated_cell_cycle_genes_regressed_out.rds")

cluster.name <- "res_1_RNA_integrated.harmony_cell_cycle_out"

# Join layers 
seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
Layers(seurat_object[["RNA"]])

# Find all markers 
all_markers <- FindAllMarkers(seurat_object, 
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
out_file <- glue("06_SILP_topDEGs/out/SILP_regressed_out_Top_100_Cluster_DE_Markers.xlsx")

# Use openxlsx::write.xlsx, which takes the named list and writes
# each element as a separate sheet (sheet name = list name, i.e., Cluster ID)
openxlsx::write.xlsx(
  x = top_markers_list,
  file = out_file,
  overwrite = TRUE # Overwrite the file if it already exists
)

################################################################################
############################## Cell cycle removed ############################## 

seurat_object <- readRDS("05_SILP_annotation/out/seurat_integrated_cell_cycle_genes_removed.rds")

cluster.name <- "res_1_RNA_integrated.harmony_cell_cycle_rm"

# Join layers 
seurat_object[["RNA"]] <- JoinLayers(seurat_object[["RNA"]])
Layers(seurat_object[["RNA"]])

# Find all markers 
all_markers <- FindAllMarkers(seurat_object, 
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
out_file <- glue("06_SILP_topDEGs/out/SILP_removed_Top_100_Cluster_DE_Markers.xlsx")

# Use openxlsx::write.xlsx, which takes the named list and writes
# each element as a separate sheet (sheet name = list name, i.e., Cluster ID)
openxlsx::write.xlsx(
  x = top_markers_list,
  file = out_file,
  overwrite = TRUE # Overwrite the file if it already exists
)

