SoupX_function <- function(data_dir, output_dir) {
  raw    <- Read10X_h5(file.path(data_dir, 'raw_feature_bc_matrix.h5') )
  filt   <- Read10X_h5(file.path(data_dir, 'filtered_feature_bc_matrix.h5') )
  srat   <- CreateSeuratObject(counts = filt)
  srat   <- SCTransform(srat, verbose = F)
  srat   <- RunPCA(srat, verbose = F)
  srat   <- RunUMAP(srat, dims = 1:30, verbose = F)
  srat   <- FindNeighbors(srat, dims = 1:30, verbose = F)
  srat   <- FindClusters(srat, verbose = F)
  meta   <- srat@meta.data
  umap   <- srat@reductions$umap@cell.embeddings
  soup   <- SoupChannel(raw, filt)
  soup   <- setClusters(soup, setNames(meta$seurat_clusters, rownames(meta)))
  soup   <- setDR(soup, umap)
  soup   <- autoEstCont(soup, doPlot = F, verbose = F)
  matrix <- adjustCounts(soup, roundToInt = T)
  write10xCounts(file.path (output_dir, 'SoupX_filtered'), matrix, overwrite = TRUE, version = '3')
}