
explore_annotation_plot <- function(seurat_obj, dimplot, marker_genes, description = "") {

  # Color UMAP by marker genes
  fp <- FeaturePlot(seurat_obj, features = marker_genes) + plot_annotation(subtitle = description)
  
  # Collapse list of marker genes to name (a single string)
  marker_genes_name <- paste(marker_genes, collapse = '_')
  
  # Save UMAP
  n_marker_genes <- length(marker_genes)
  
  if (n_marker_genes == 1){
    
    ggsave(glue("05_annotation/plot/UMAP_{marker_genes_name}.pdf"), plot = fp, width = 7, height = 6)
    
  } else if (n_marker_genes == 2){
    
    ggsave(glue("05_annotation/plot/UMAP_{marker_genes_name}.pdf"), plot = fp, width = 12, height = 6)
    
  } else if (n_marker_genes == 3){
    
    ggsave(glue("05_annotation/plot/UMAP_{marker_genes_name}.pdf"), plot = fp, width = 12, height = 12)
    
  } else {
    
    ggsave(glue("05_annotation/plot/UMAP_{marker_genes_name}.pdf"), plot = fp)
    
  }
  
  
  # Wrap fp with UMAP colored by clusters and save 
  if (n_marker_genes == 1){
    
    wp <- dimplot + fp + plot_annotation(subtitle = description)
    ggsave(glue("05_annotation/plot/UMAP_clusters_{marker_genes_name}.pdf"), plot = wp, width = 12, height = 6)
    
  } else if (n_marker_genes == 2){
    
    wp <- (dimplot + ggplot()) / fp + plot_annotation(subtitle = description)
    ggsave(glue("05_annotation/plot/UMAP_clusters_{marker_genes_name}.pdf"), plot = wp, width = 14, height = 12)
    
  } else if (n_marker_genes == 3){
    
    wp <- dimplot + fp + plot_annotation(subtitle = description)
    ggsave(glue("05_annotation/plot/UMAP_clusters_{marker_genes_name}.pdf"), plot = wp, width = 14, height = 7)
    
  } else {
    
    wp <- dimplot + fp + plot_annotation(subtitle = description)
    ggsave(glue("05_annotation/plot/UMAP_clusters_{marker_genes_name}.pdf"), plot = wp) 
    
  }
  
  
}

