
library(zellkonverter)
adata <- readH5AD("~/Downloads/Myeloid_allgenes.h5ad")
adata

assayNames(adata)


x <- assay(adata, "X")
all(x == round(x)) # check if counts 
