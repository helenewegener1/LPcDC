# Load libraries 
library(SeuratObject)
library(Seurat)
library(dplyr)
library(stringr)
library(glue)
library(ggplot2)

# Load data
seurat_object <- readRDS("34_STROMAL_integrate/out/STROMAL_seurat_integrated_v5_RNA.rds")

res <- c(1, 1.5)
method <- c("harmony", "cca")
# method <- "cca"

# Create all combinations
combos <- expand.grid(method = method, res = res)

# Generate cluster names
cluster.names <- glue("RNA_{combos$method}_clusters_res.{combos$res}") %>% as.character()

for (cluster.name in cluster.names){
  
  # cluster.name <- "RNA_cca_clusters_res.1"
  
  # Join layers before DE analysis 
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
  out_file <- glue("36_STROMAL_topDEGs/out/STROMAL_Top_100_DE_Markers_{cluster.name}.xlsx")
  
  # Use openxlsx::write.xlsx, which takes the named list and writes
  # each element as a separate sheet (sheet name = list name, i.e., Cluster ID)
  openxlsx::write.xlsx(
    x = top_markers_list,
    file = out_file,
    overwrite = TRUE # Overwrite the file if it already exists
  )
  
}

