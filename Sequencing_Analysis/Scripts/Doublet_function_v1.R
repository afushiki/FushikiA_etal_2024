doublet_exclusion <- function(x, doublets_rate) {
  x    <- CreateSeuratObject(counts = x, min.cells = 3, min.features = 200)
  x    <- PercentageFeatureSet(x, pattern = '^(MT|Mt|mt)-', col.name = 'percent_mt')
  x    <- SCTransform (x, vars.to.regress = 'percent_mt', verbose = FALSE)
  x    <- RunPCA(x, verbose = F, npcs = 20)
  x    <- RunUMAP(x, dims = 1:10, metric = 'correlation', verbose = F)
  nExp <- round(ncol(x) * doublets_rate) 
  x    <- doubletFinder(x, pN = 0.25, pK = 0.09, sct=TRUE, nExp = nExp, PCs = 1:10)
  df   <- colnames(x@meta.data)[grepl('DF.classification', colnames(x@meta.data))]
  x    <- x[, x@meta.data[, df] == 'Singlet']
}