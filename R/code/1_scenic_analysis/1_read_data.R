suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratDisk)
  library(SeuratObject)
  library(SeuratData)
  library(loomR)
  library(hdf5r)
  library(LoomExperiment)
  library(SCopeLoomR)
  library(patchwork)
  library(dplyr)
  library(purrr)
  library(graphics)
  library(grDevices)
  library(utils)
  library(datasets)
  library(methods)
  library(base)
  library(stats)
  library(tidyverse)
  library(ggplot2)
})
rm(list = ls())

# Convert the data file to H5Seurat format
Convert(
  'adata_body_S_v1.0.h5ad',
  dest = "h5seurat",
  overwrite = TRUE
)

### Seurat Import
afca_headBody_seurat <- LoadH5Seurat("adata_headBody_S_v1.0.h5seurat",meta.data = T)
afca_head_seurat <- LoadH5Seurat("adata_head_S_v1.0.h5seurat",meta.data = T)
afca_body_seurat <- LoadH5Seurat("adata_body_S_v1.0.h5seurat",meta.data = T)

# Save the Seurat objects to RDS files
saveRDS(afca_headBody_seurat, "afca_headBody_seurat.RDS")
saveRDS(afca_head_seurat, "afca_head_seurat.RDS")
saveRDS(afca_body_seurat, "afca_body_seurat.RDS")
