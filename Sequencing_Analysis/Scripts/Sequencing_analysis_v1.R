#---Clear---#
rm(list = ls())
cat('\014')
gc() 

suppressPackageStartupMessages({
  library(Seurat)
  library(SeuratWrappers)
  library(scCustomize)
  library(tidyverse)
  library(cowplot) 
  library(SoupX)
  library(DoubletFinder)
  library(reticulate)
  library(patchwork)
  library(SingleCellExperiment)
  library(miloR)
  library(lme4)
  library(anndata)
  library(gprofiler2)
  library(gt)
  library(gtsummary) 
  library(DropletUtils)
  library(knitr)
  library(RColorBrewer)
  library(ggpubr)
  library(rstatix)
  library(lmerTest)  
  library(emmeans)  
  library(readr)
  library(ggplot2)
})

options(future.globals.maxSize = 100000 * 1024^2)

source('./Scripts/Seurat/SoupX_function_v1.R')
source('./Scripts/Seurat/Doublet_function_v1.R')
source('./Scripts/Seurat/MASC_analysis_v1.R')
source('./Scripts/Seurat/Milo_analysis_v1.R')
source('./Scripts/Seurat/GO_analysis_v1.R')
source('./Scripts/Seurat/Proportion_analysis_df_stat_v1.R')
source('./Scripts/Seurat/Proportion_analysis_plot_v1.R')

filename           <- 'Seq_PD_model_mice'
data_dir           <- './RNAseq_data'
results            <- file.path(data_dir, paste0('RNAseq_analysis', filename))
output             <- file.path(results, 'output')
go_analysis_output <- file.path(output, 'GO_analysis_output')
dir.create(results)
dir.create(file.path(results, 'output'))
dir.create(file.path(output, 'GO_analysis_output'))

cell_id_list <- c('RA004', 'RA012', 'RA017', 'RA018', 'RA005', 'RA006', 'RA013', 'RA021', 'RA015', 'RA016', 'RA019', 'RA020', 
                  'RA034', 'RA047', 'RA039', 'RA049', 'RA043', 'RA048', 'RA040', 'RA041', 'RA035', 'RA036', 'RA044', 'RA050',
                  'RA026', 'RA027', 'RA030', 'RA031', 'RA025', 'RA032', 'RA046', 'RA051', 'RA022', 'RA023', 'RA028', 'RA029')

cr_list_data <- list()
for(i in 1:length(cell_id_list)){
  cr_list_data <- append(cr_list_data, file.path(data_dir, paste0(cell_id_list[i], '_cellranger_count_outs')) )
}

cr_list_results <- list()
for(i in 1:length(cell_id_list)){
  cr_list_results <- append(cr_list_results, file.path(data_dir, paste0(cell_id_list[i], '_cellranger_count_outs')) )
  dir.create(cr_list_results[[i]])
}

# QC
mapply(SoupX_function, cr_list_data, cr_list_results)
soupx_path <- lapply(cr_list_results, function(x) file.path (x, 'SoupX_filtered'))
soupx_list <- lapply(soupx_path, Read10X)

doublets_rate_list     <- c(0.023, 0.069, 0.084, 0.061, 0.069, 0.061, 0.135, 0.046, 0.031, 0.046, 0.061, 0.054, 
                            0.054, 0.054, 0.061, 0.054, 0.054, 0.054, 0.054, 0.046, 0.084, 0.046, 0.084, 0.054,
                            0.061, 0.069, 0.054, 0.076, 0.039, 0.069, 0.054, 0.046, 0.054, 0.069, 0.084, 0.054) 
doublet_excluded_list  <- mapply(function(X,Y) doublet_exclusion(X, Y), soupx_list, doublets_rate_list)

seu_list <- lapply(doublet_excluded_list, function(x) CreateSeuratObject(x@assays$RNA$counts, min.cells = 3, min.features=200))
seu_list <- lapply(1:length(seu_list), function(x) {
  seu_list[[x]] <- PercentageFeatureSet(seu_list[[x]], pattern = '^(MT|Mt|mt)-', col.name = 'percent_mt')
  seu_list[[x]] <- PercentageFeatureSet(seu_list[[x]], pattern = '^(RP[LS]|rp[ls])', col.name = 'percent_ribo')
  seu_list[[x]] <- PercentageFeatureSet(seu_list[[x]], pattern = '^(HB[^(P|E|S)]|Hb[^(p|e|s)])', col.name = 'percent_hemo')
})
seu_list <- lapply(seu_list, function(x) {x <- subset(x, subset = nFeature_RNA > 1000 & nFeature_RNA < 7500 & nCount_RNA < 40000 & percent_mt < 1 & percent_ribo < 1 & percent_hemo < 1)})

# For integration
seu_list <- lapply(X = seu_list, FUN = function(x) {
  x <- NormalizeData(x, verbose = FALSE)
  x <- FindVariableFeatures(x, verbose = FALSE)
})
features <- SelectIntegrationFeatures(object.list = seu_list)
seu_list <- lapply(X = seu_list, FUN = function(x) {
  x <- ScaleData(x, features = features, verbose = FALSE)
  x <- RunPCA(x, features = features, verbose = FALSE)
})
IntData <- Merge_Seurat_List(seu_list, project = 'IntData', add.cell.ids = cell_id_list)
IntData <- AddMetaData(IntData, metadata = substr(rownames(IntData@meta.data), 1, 5), col.name = 'cell_ids')

# For annotation
IntData[['RNA']] <- JoinLayers(IntData[['RNA']]) 
dataIn <- LayerData(object = IntData, assay = 'RNA', layer = 'counts') 
dataQC <- dataIn[,colSums(dataIn)>250] 
dataQCt = Matrix::t(dataQC)
ad <- AnnData(X = dataQCt,
              obs = data.frame(group = rownames(dataQCt), row.names = rownames(dataQCt)),
              var = data.frame(type = colnames(dataQCt), row.names = colnames(dataQCt)))
write_h5ad(ad,file.path(results, paste0(filename, '.h5ad')),compression='gzip')

# Go to the website of Mapmycell and run h5ad files

mapping <- read.csv(file.path(results, 'annotation/mapmycells-output-af_1715006792868.csv'),comment.char='#')
mapping$class_new <- mapping$class_name
mapping$class_new[!is.element(mapping$class_name,names(head(-sort(-table(mapping$class_name)),15)))] = 'other' # to visualize microglia group, increaased from 8 to 15
# mapping$class_new[!is.element(mapping$class_name,names(head(-sort(-table(mapping$class_name)),8)))] = 'other' # original
mapping$subclass_new <- mapping$subclass_name
mapping$subclass_new[!is.element(mapping$subclass_name,names(head(-sort(-table(mapping$subclass_name)),30)))] = 'other'  # to visualize microglia group, increaased from 20 to 30
# mapping$subclass_new[!is.element(mapping$subclass_name,names(head(-sort(-table(mapping$subclass_name)),20)))] = 'other' # original
rownames(mapping) <- mapping$cell_id
mapping <- mapping[colnames(dataQC),]
dataSeurat <- CreateSeuratObject(counts = dataQC, meta.data = mapping)
dataSeurat <- AddMetaData(dataSeurat, metadata = substr(rownames(dataSeurat@meta.data), 1, 5), col.name = 'cell_ids')

# Add groups
MP_8wk     <- subset(dataSeurat, subset = cell_ids %in% c('RA004', 'RA012', 'RA017', 'RA018')); MP_8wk[['groups']]     <- 'a_MP_8wk'
MP_16wk    <- subset(dataSeurat, subset = cell_ids %in% c('RA005', 'RA006', 'RA013', 'RA021')); MP_16wk[['groups']]    <- 'b_MP_16wk'
MP_24wk    <- subset(dataSeurat, subset = cell_ids %in% c('RA015', 'RA016', 'RA019', 'RA020')); MP_24wk[['groups']]    <- 'c_MP_24wk'
MP_lm_8wk  <- subset(dataSeurat, subset = cell_ids %in% c('RA034', 'RA047', 'RA039', 'RA049')); MP_lm_8wk[['groups']]  <- 'd_MP_lm_8wk'
MP_lm_16wk <- subset(dataSeurat, subset = cell_ids %in% c('RA043', 'RA048', 'RA040', 'RA041')); MP_lm_16wk[['groups']] <- 'e_MP_lm_16wk'
MP_lm_24wk <- subset(dataSeurat, subset = cell_ids %in% c('RA035', 'RA036', 'RA044', 'RA050')); MP_lm_24wk[['groups']] <- 'f_MP_lm_24wk'
WT_8wk     <- subset(dataSeurat, subset = cell_ids %in% c('RA026', 'RA027', 'RA030', 'RA031')); WT_8wk[['groups']]     <- 'g_WT_8wk'
WT_16wk    <- subset(dataSeurat, subset = cell_ids %in% c('RA025', 'RA032', 'RA046', 'RA051')); WT_16wk[['groups']]    <- 'h_WT_16wk'
WT_24wk    <- subset(dataSeurat, subset = cell_ids %in% c('RA022', 'RA023', 'RA028', 'RA029')); WT_24wk[['groups']]    <- 'i_WT_24wk'

IntData    <- merge(MP_8wk, y = c(MP_16wk, MP_24wk, MP_lm_8wk, MP_lm_16wk, MP_lm_24wk, WT_8wk, WT_16wk, WT_24wk), project = 'IntData')
DefaultAssay(IntData) <- 'RNA'

# Basic process 
IntData <- PercentageFeatureSet(IntData, pattern = '^(MT|Mt|mt)-', col.name = 'percent_mt')
IntData <- NormalizeData(IntData, verbose = FALSE)
IntData <- FindVariableFeatures(IntData, verbose = FALSE)
IntData <- ScaleData(IntData, vars.to.regress = 'percent_mt', verbose = FALSE)
IntData <- RunPCA(IntData, verbose = FALSE)
# Unintegrated
IntData <- FindNeighbors(IntData, reduction = 'pca', dims = 1:30)
IntData <- FindClusters(IntData, resolution = 0.8, cluster.name = 'unintegrated_clusters')
IntData <- RunUMAP(IntData, reduction = 'pca', dims = 1:30, reduction.name = 'umap.unintegrated')
# RPCA integration
IntData <- IntegrateLayers(IntData, method = RPCAIntegration, orig.reduction = 'pca', new.reduction = 'integrated.rpca', verbose = FALSE)
IntData <- FindNeighbors(IntData, reduction = 'integrated.rpca', dims = 1:30)
IntData <- FindClusters(IntData, resolution = 0.8, cluster.name = 'rpca_clusters')
IntData <- RunUMAP(IntData, reduction = 'integrated.rpca', dims = 1:30, reduction.name = 'umap.rpca')
# Harmony integration
IntData <- IntegrateLayers(IntData, method = HarmonyIntegration, orig.reduction = 'pca', new.reduction = 'harmony', verbose = FALSE)
IntData <- FindNeighbors(IntData, reduction = 'harmony', dims = 1:30)
IntData <- FindClusters(IntData, resolution = 0.8, cluster.name = 'harmony_clusters')
IntData <- RunUMAP(IntData, reduction = 'harmony', dims = 1:30, reduction.name = 'umap.harmony')
# scVI integration
IntData <- IntegrateLayers(IntData, method = scVIIntegration, new.reduction = 'integrated.scvi', conda_env = 'C:/Users/anaconda3/envs/scvi-gpu', verbose = FALSE)
IntData <- FindNeighbors(IntData, reduction = 'integrated.scvi', dims = 1:30)
IntData <- FindClusters(IntData, resolution = 0.8, cluster.name = 'scvi_clusters')
IntData <- RunUMAP(IntData, reduction = 'integrated.scvi', dims = 1:30, reduction.name = 'umap.scvi')

IntData[['RNA']] <- JoinLayers(IntData[['RNA']]) 
IntData@reductions$umap <- IntData@reductions$umap.harmony # this is for milo analysis
# The original annotation contains spaces or hyphens, so please replace them with underscores.
IntData@meta.data$class_new      <- gsub(' ', '_', IntData@meta.data$class_new) 
IntData@meta.data$class_new      <- gsub('-', '_', IntData@meta.data$class_new) 
IntData@meta.data$subclass_new   <- gsub(' ', '_', IntData@meta.data$subclass_new) 
IntData@meta.data$subclass_new   <- gsub('-', '_', IntData@meta.data$subclass_new )
IntData@meta.data$supertype_name <- gsub(' ', '_', IntData@meta.data$supertype_name) 
IntData@meta.data$supertype_name <- gsub('-', '_', IntData@meta.data$supertype_name)
IntData@meta.data$cluster_name   <- gsub(' ', '_', IntData@meta.data$cluster_name) 
IntData@meta.data$cluster_name   <- gsub('-', '_', IntData@meta.data$cluster_name)
IntData <- AddMetaData(IntData, metadata = IntData@meta.data$class_new, col.name = 'class_new') # it's for MASC analysis
IntData <- AddMetaData(IntData, metadata = IntData@meta.data$subclass_new, col.name = 'subclass_new') 
IntData <- AddMetaData(IntData, metadata = IntData@meta.data$supertype_name, col.name = 'supertype_name') 
IntData <- AddMetaData(IntData, metadata = IntData@meta.data$cluster_name, col.name = 'cluster_name') 

# Color scheme
color_pal_int1 <- DiscretePalette_scCustomize(num_colors = n_distinct(IntData@meta.data$class_new), palette = 'polychrome', shuffle_pal = TRUE)
color_pal_int2 <- DiscretePalette_scCustomize(num_colors = n_distinct(IntData@meta.data$subclass_new), palette = 'varibow', shuffle_pal = TRUE)

# Create images
png(file.path(output, 'IntData_MitoPark_harmony_cluster_class_subclass_w_legend.png'), width = 15, height = 5, units = 'in', res = 300)
p1 <- DimPlot(IntData, reduction = 'umap.harmony', group.by='class_new', cols=color_pal_int1, label=FALSE, label.size = 4, repel=TRUE, raster=FALSE)  + NoAxes()
p2 <- DimPlot(IntData, reduction = 'umap.harmony', group.by='subclass_new', cols=color_pal_int2, label=FALSE, label.size = 3, repel=TRUE, raster=FALSE) +
   theme(legend.text=element_text(size=8), legend.key.size = unit(0.5, 'line')) + guides(color=guide_legend(ncol =1, override.aes = list(size = 4))) + NoAxes()
p1 + p2
dev.off()

png(file.path(output, 'IntData_Th_MitoPark_Cluster_Th_Slc17a6_Slc32a1.png'), width = 18, height = 4, units = 'in', res = 300)
p1 <- DimPlot(IntData, reduction = 'umap.harmony', group.by='class_new', cols=color_pal_int1, label=TRUE, label.size = 4, repel=TRUE) + NoAxes()
p2 <- FeaturePlot(IntData, features = 'Th', reduction = 'umap.harmony', max.cutoff = 'q95', cols = c('lightgrey', 'red'), raster=FALSE) + 
  theme(plot.title = element_text(size = 15, face = 'bold'),legend.position = c(0.85, 0.2), legend.text=element_text(size=12), legend.key.size = unit(1, 'line')) +
  scale_color_gradient(low = 'lightgrey', high = 'red',  breaks = c(0, 1, 2)) + NoAxes()
p3 <- FeaturePlot(IntData, features = 'Slc17a6', reduction = 'umap.harmony', max.cutoff = 'q95', cols = c('lightgrey', 'red'), raster=FALSE) + 
  theme(plot.title = element_text(size = 15, face = 'bold'),legend.position = c(0.85, 0.2), legend.text=element_text(size=12), legend.key.size = unit(1, 'line')) +
  scale_color_gradient(low = 'lightgrey', high = 'red',  breaks = c(0, 1, 2)) + NoAxes()
p4 <- FeaturePlot(IntData, features = 'Slc32a1', reduction = 'umap.harmony', max.cutoff = 'q95', cols = c('lightgrey', 'red'), raster=FALSE) + 
  theme(plot.title = element_text(size = 15, face = 'bold'),legend.position = c(0.85, 0.2), legend.text=element_text(size=12), legend.key.size = unit(1, 'line')) +
  scale_color_gradient(low = 'lightgrey', high = 'red',  breaks = c(0, 1, 2)) + NoAxes()
p1|p2|p3|p4
dev.off()

png(file.path(output, 'IntData_MitoPark_split_by_group.png'), width = 15, height = 12, units = 'in', res = 300)
DimPlot(IntData, reduction = 'umap.harmony', group.by='class_new', split.by = 'groups', ncol = 3, cols=color_pal_int1, label=FALSE, label.size = 4, repel=TRUE) +
  theme(legend.position = c(1, 0.5), legend.text=element_text(size=20), legend.key.size = unit(1.5, 'line')) +
  theme(plot.margin=unit(c(0,3,0,0), 'in')) + NoAxes()
dev.off()

png(file.path(output, 'IntData_MitoPark_cluster_genes_expression.png'), width = 12, height = 8, units = 'in', res = 300)
FeaturePlot(IntData, features = c('Rbfox3', 'Mog', 'Pdgfra', 'Ptprc', 'Gja1', 'Nostrin'),
            max.cutoff = 'q95', cols = c('lightgrey', 'red'), order=T, raster=FALSE, reduction = 'umap.harmony',ncol=3) & NoAxes() & theme(legend.position = c(0.85, 0.2), legend.text=element_text(size=7), legend.key.size = unit(0.5, 'line')) 
dev.off()

# Subsetting by dopamine cluster with Allen's annotation
IntData_Th <- subset(IntData, class_new == '21_MB_Dopa')
DefaultAssay(IntData_Th) <- 'RNA'
res_cluster = 0.8
IntData_Th <- ScaleData(IntData_Th, vars.to.regress = 'percent_mt', verbose = FALSE) 
IntData_Th <- RunPCA(IntData_Th, verbose = FALSE)
IntData_Th <- FindNeighbors(IntData_Th, dims = 1:30)
IntData_Th <- FindClusters(IntData_Th, resolution = res_cluster, cluster.name = 'umap_clusters')
IntData_Th <- RunUMAP(IntData_Th, dims = 1:30, reduction.name = 'umap')

# Rename to short annotation
IntData_Th@meta.data$supertype_name <- gsub('_SNc_VTA_RAmb_Foxa1_Dopa_([0-9])', '', IntData_Th@meta.data$supertype_name) 
IntData_Th@meta.data$cluster_name   <- gsub('_SNc_VTA_RAmb_Foxa1_Dopa_([0-9])', '', IntData_Th@meta.data$cluster_name) 
IntData_Th <- AddMetaData(IntData_Th, metadata = IntData_Th@meta.data$supertype_name, col.name = 'supertype_name_short') 
IntData_Th <- AddMetaData(IntData_Th, metadata = IntData_Th@meta.data$cluster_name, col.name = 'cluster_name_short') 
 
# For downstream analysis
IntData_Th_supertype <- IntData_Th; Idents(IntData_Th_supertype) <- IntData_Th_supertype@meta.data$supertype_name_short
IntData_Th_cluster   <- IntData_Th; Idents(IntData_Th_cluster)   <- IntData_Th_cluster@meta.data$cluster_name_short
color_pal1 <- DiscretePalette_scCustomize(num_colors = length(levels(Idents(IntData_Th_supertype))), palette = 'polychrome')
color_pal2 <- DiscretePalette_scCustomize(num_colors = length(levels(Idents(IntData_Th_cluster))), palette = 'varibow')

# Marker analysis
IntData_Th_supertype_markers <- FindAllMarkers(IntData_Th_supertype, min.pct = 0.25, logfc.threshold = 0.25, assay='RNA', verbose = F)
IntData_Th_cluster_markers   <- FindAllMarkers(IntData_Th_cluster, min.pct = 0.25, logfc.threshold = 0.25, assay='RNA', verbose = F)
top5_IntData_Th_supertype    <- Extract_Top_Markers(IntData_Th_supertype_markers, num_genes = 5, rank_by = 'avg_log2FC', named_vector = FALSE, make_unique = TRUE)
top3_IntData_Th_supertype    <- Extract_Top_Markers(IntData_Th_supertype_markers, num_genes = 3, rank_by = 'avg_log2FC', named_vector = FALSE, make_unique = TRUE)
top1_IntData_Th_cluster      <- Extract_Top_Markers(IntData_Th_cluster_markers, num_genes = 1, rank_by = 'avg_log2FC', named_vector = FALSE, make_unique = TRUE)

# Create images
png(file.path(output, paste0('IntData_Th_supertype_markers_heatmap_top5_annotated.png')), width = 18, height =16, units = 'in', res = 300)
p <- DoHeatmap(IntData_Th_supertype, top5_IntData_Th_supertype, size = 10) + theme(text = element_text(size = 30), legend.text=element_text(size=25), legend.key.size = unit(3, 'line')) + 
  guides(colour=FALSE) # top5_IntData_Th_supertype$gene
print(p)
dev.off()

png(file.path(output, 'IntData_MitoPark_Th_subcluster_genes_expression.png'), width = 16, height = 9, units = 'in', res = 300)
FeaturePlot(IntData_Th, features = c('Th','Ddc', 'Sox6', 'Aldh1a1','Anxa1', 'Calb1', 'Slc17a6','Otx2'),
            max.cutoff = 'q95', cols = c('lightgrey', 'red'), order=T, raster=FALSE, reduction = 'umap',ncol=4) &
  theme(plot.title = element_text(size = 25, face = 'bold'), legend.position = c(0.85, 0.1), legend.text=element_text(size=10), legend.key.size = unit(1.0, 'line'), 
        strip.text.x = element_text(size = 20, face = 'bold')) & NoAxes() 
dev.off()

png(file.path(output, 'IntData_MitoPark_subset_allen_taxonomy_all.png'), width = 10, height = 10, units = 'in', res = 300)
p1 <- DimPlot(IntData_Th, reduction = 'umap', group.by='class_name', cols=color_pal1, label=TRUE, repel = TRUE, label.size = 10) + 
  labs(title = 'Class') + theme(plot.title = element_text(hjust = 0.5)) + NoLegend() + NoAxes()
p2 <- DimPlot(IntData_Th, reduction = 'umap', group.by='subclass_name', cols=color_pal2, label=TRUE, repel = TRUE, label.size = 6)  + 
  labs(title = 'Subclass') + theme(plot.title = element_text(hjust = 0.5)) + NoLegend() + NoAxes()
p3 <- DimPlot(IntData_Th, reduction = 'umap', group.by='supertype_name', cols=color_pal1, label=TRUE, repel = TRUE, label.size = 6) + 
  labs(title = 'Supertype') + theme(plot.title = element_text(hjust = 0.5)) + NoLegend() + NoAxes()
p4 <- DimPlot(IntData_Th, reduction = 'umap', group.by='cluster_name', cols=color_pal2, label=TRUE, repel = TRUE, label.size = 4)  + 
  labs(title = 'Cluster') + theme(plot.title = element_text(hjust = 0.5)) + NoLegend() + NoAxes()
(p1|p2)/(p3|p4)
dev.off()

# Re-annotate Th cluster
supertype_markers        <- c('Megf11', 'Etv1', 'Plekhg1', 'Sox6', 'Sema5b', 'Cck', 'Synpr', 'Lypd1')
cluster_markers          <- top1_IntData_Th_cluster
new_ids_supertype        <- paste0(levels(Idents(IntData_Th_supertype)), '_', supertype_markers) 
names(new_ids_supertype) <- levels(IntData_Th_supertype)
IntData_Th_supertype     <- RenameIdents(IntData_Th_supertype, new_ids_supertype)
IntData_Th_supertype@meta.data$supertype_name <- paste(Idents(IntData_Th_supertype)) 
new_ids_cluster          <- paste0(levels(Idents(IntData_Th_cluster)), '_', cluster_markers) 
names(new_ids_cluster)   <- levels(IntData_Th_cluster)
IntData_Th_cluster       <- RenameIdents(IntData_Th_cluster, new_ids_cluster)
IntData_Th_cluster@meta.data$cluster_name <- paste(Idents(IntData_Th_cluster)) 

# Create images
png(file.path(output, paste0('IntData_Th_supertype_MitoPark_marker_genes_expression.png')), width = 12, height = 6, units = 'in', res = 300)
FeaturePlot(IntData_Th_supertype, features = supertype_markers,  max.cutoff = 'q95', cols = c('lightgrey', 'red'), order=T, raster=FALSE, reduction = 'umap', ncol=4) & NoAxes() & 
  theme(plot.title = element_text(size = 15, face = 'bold'),legend.position = c(0.9, 0.2), legend.text=element_text(size=12), legend.key.size = unit(0.7, 'line'))  &
  scale_color_gradient(low = 'grey', high = 'red',  breaks = c(0, 1, 2))  
dev.off()

png(file.path(output, paste0('IntData_Th_supertype_MitoPark_Clustered_DotPlot.png')), width = 8, height = 8, units = 'in', res = 300)
p <- Clustered_DotPlot(IntData_Th_supertype, features = top3_IntData_Th_supertype, column_label_size = 15, row_label_size = 15, x_lab_rotate=45, k = length(levels(Idents(IntData_Th_supertype)))) # k = number of cluster
print(p)
dev.off()

png(file.path(output, 'IntData_MitoPark_subset_supertype_marker_annotated.png'), width = 5, height = 6, units = 'in', res = 300)
p <- DimPlot(IntData_Th_supertype, reduction = 'umap', group.by='supertype_name', cols=color_pal1, label=TRUE, repel = TRUE, label.size=6) + 
  theme(legend.text=element_text(size=9), legend.key.size = unit(1, 'line'), legend.position='bottom') + 
  guides(color=guide_legend(ncol =4, override.aes = list(size = 6))) + NoAxes()
p
dev.off()

png(file.path(output, paste0('IntData_Th_supertype_MitoPark_Stacked_top1_VlnPlot.png')), width = 5, height = 5, units = 'in', res = 300)
Stacked_VlnPlot(IntData_Th_supertype, supertype_markers, x_lab_rotate = TRUE, 
                colors_use = DiscretePalette_scCustomize(num_colors = length(levels(Idents(IntData_Th_supertype))), palette = 'polychrome')) +
  ggtitle('Identity on x-axis') + xlab('Identity') + ylab('Expression Level')
dev.off()

png(file.path(output, 'IntData_MitoPark_Th_supertype_cluster_reannotated.png'), width = 5, height = 5.5, units = 'in', res = 300)
p <- DimPlot(IntData_Th, reduction = 'umap', group.by='supertype_name', split.by = 'groups', cols=color_pal1, ncol=3) +
  theme(strip.text.x = element_text(size = 10, face = 'bold')) + ggtitle(NULL) + NoLegend() + NoAxes()
p
dev.off()

# Differential abundance tests: Milo and MASC
IntData_Th_supertype <- AddMetaData(IntData_Th_supertype, metadata = substr(rownames(IntData_Th_supertype@meta.data), 1, 5), col.name = 'cell_ids')
IntData_Th_supertype <- AddMetaData(IntData_Th_supertype, metadata = IntData_Th_supertype@active.ident, col.name = 'marker') 
IntData_Th_cluster   <- AddMetaData(IntData_Th_cluster, metadata = substr(rownames(IntData_Th_cluster@meta.data), 1, 5), col.name = 'cell_ids')
IntData_Th_cluster   <- AddMetaData(IntData_Th_cluster, metadata = IntData_Th_cluster@active.ident, col.name = 'marker') 
IntData_Th_supertype_list <- SplitObject(IntData_Th_supertype, split.by = 'groups')
IntData_Th_cluster_list   <- SplitObject(IntData_Th_cluster, split.by = 'groups')

MP_16wk_w_lm_Th_supertype  <- merge(IntData_Th_supertype_list[[5]], y = IntData_Th_supertype_list[[2]], project = 'MP_16wk_w_lm_Th_supertype',  merge.dr = TRUE)
MP_24wk_w_lm_Th_supertype  <- merge(IntData_Th_supertype_list[[6]], y = IntData_Th_supertype_list[[3]], project = 'MP_24wk_w_lm_Th_supertype',  merge.dr = TRUE)
MP_16wk_w_8wk_Th_supertype <- merge(IntData_Th_supertype_list[[2]], y = IntData_Th_supertype_list[[1]], project = 'MP_16wk_w_8wk_Th_supertype', merge.dr = TRUE)
MP_16wk_w_lm_Th_supertype@meta.data$groups <- gsub('e_MP_lm_16wk', '1_MP_lm_16wk', MP_16wk_w_lm_Th_supertype@meta.data$group)
MP_16wk_w_lm_Th_supertype@meta.data$groups <- gsub('b_MP_16wk', '2_MP_16wk', MP_16wk_w_lm_Th_supertype@meta.data$group)
MP_24wk_w_lm_Th_supertype@meta.data$groups <- gsub('f_MP_lm_24wk', '3_MP_lm_24wk', MP_24wk_w_lm_Th_supertype@meta.data$group)
MP_24wk_w_lm_Th_supertype@meta.data$groups <- gsub('c_MP_24wk', '4_MP_24wk', MP_24wk_w_lm_Th_supertype@meta.data$group)

MP_16wk_w_lm_Th_cluster  <- merge(IntData_Th_cluster_list[[5]], y = IntData_Th_cluster_list[[2]], project = 'MP_16wk_w_lm_Th_cluster',  merge.dr = TRUE)
MP_24wk_w_lm_Th_cluster  <- merge(IntData_Th_cluster_list[[6]], y = IntData_Th_cluster_list[[3]], project = 'MP_24wk_w_lm_Th_cluster',  merge.dr = TRUE)
MP_16wk_w_8wk_Th_cluster <- merge(IntData_Th_cluster_list[[2]], y = IntData_Th_cluster_list[[1]], project = 'MP_16wk_w_8wk_Th_cluster',  merge.dr = TRUE)
MP_16wk_w_lm_Th_cluster@meta.data$groups <- gsub('e_MP_lm_16wk', '1_MP_lm_16wk', MP_16wk_w_lm_Th_cluster@meta.data$group)
MP_16wk_w_lm_Th_cluster@meta.data$groups <- gsub('b_MP_16wk', '2_MP_16wk', MP_16wk_w_lm_Th_cluster@meta.data$group)
MP_24wk_w_lm_Th_cluster@meta.data$groups <- gsub('f_MP_lm_24wk', '3_MP_lm_24wk', MP_24wk_w_lm_Th_cluster@meta.data$group)
MP_24wk_w_lm_Th_cluster@meta.data$groups <- gsub('c_MP_24wk', '4_MP_24wk', MP_24wk_w_lm_Th_cluster@meta.data$group)

# Milo analysis: supertype
milo_analysis(MP_16wk_w_lm_Th_supertype, 'MP_16wk_w_lm_Th_supertype', 'umap', 'supertype_name', 0.7, color_pal1, 2, TRUE, 8, c(-8, 8), c(-8, 8), c(-8, 0, 8), 28)
milo_analysis(MP_24wk_w_lm_Th_supertype, 'MP_24wk_w_lm_Th_supertype', 'umap', 'supertype_name', 0.7, color_pal1, 2, TRUE, 8, c(-8, 8), c(-8, 8), c(-8, 0, 8), 28)
milo_analysis(MP_16wk_w_8wk_Th_supertype, 'MP_16wk_w_8wk_Th_supertype', 'umap', 'supertype_name', 0.7, color_pal1, 2, TRUE, 8, c(-3, 3), c(-3, 3), c(-3, 0, 3), 28)
milo_analysis(MP_16wk_w_lm_Th_cluster, 'MP_16wk_w_lm_Th_cluster', 'umap', 'cluster_name', 0.2, color_pal2, 2, TRUE, 4,  c(-8, 8), c(-8, 8), c(-8, 0, 8),xaxis_size=18)
milo_analysis(MP_24wk_w_lm_Th_cluster, 'MP_24wk_w_lm_Th_cluster', 'umap', 'cluster_name', 0.2, color_pal2, 2, TRUE, 4,  c(-8, 8), c(-8, 8), c(-8, 0, 8),xaxis_size=18)
milo_analysis(MP_16wk_w_8wk_Th_cluster , 'MP_16wk_w_8wk_Th_cluster', 'umap', 'cluster_name', 0.2, color_pal2, 2, TRUE, 4, c(-3, 3), c(-3, 3), c(-3, 0, 3),xaxis_size=18)

# MASC analysis
MASC_analysis(IntData, 'IntData', 'class_new', xlim=2, xlim_log=2, angle=45, width_bubble = 20, height_bubble = 12, element_text_y=12) # xlim for dotplots
MASC_analysis(IntData, 'IntData', 'subclass_new', xlim=2, xlim_log=3, angle=45, width_bubble = 30, height_bubble = 22, element_text_y=10) # xlim for dotplots
MASC_analysis(IntData_Th_supertype, 'IntData_Th_supertype', 'marker', xlim=4, xlim_log=3, angle=45, width_bubble = 15, height_bubble = 10, element_text_y=12) # xlim for dotplots
MASC_analysis(IntData_Th_cluster, 'IntData_Th_cluster', 'marker', xlim=5, xlim_log=5, angle=45, width_bubble = 15, height_bubble = 32, element_text_y=10) # xlim for dotplots

##--Proportional plot update---###
my_comparisons            <- list(c('MitoPark 8wks', 'MitoPark 16wks', 'MitoPark 8wks', 'MitoPark 8wks', 'MitoPark 16wks', 'MitoPark 24wks'), 
                                  c('MitoPark 16wks', 'MitoPark 24wks', 'MitoPark 24wks', 'Littermate 8wks', 'Littermate 16wks', 'Littermate 24wks'))

IntData_class_prop          <- Prorpotiona_analysis_create_df(IntData, 'IntData', 'class_new') 
IntData_class_prop_filtered <- IntData_class_prop$Prop_noB6 %>% filter(!str_starts(cluster, c('01_IT_|04_DG_|12_HY_|23_P_|26_P_|30_Astro_|34_Immune|other')))   
stats_list                  <- Prorpotiona_analysis_stat(IntData_class_prop_filtered, 'class_new', my_comparisons, 'Prop_noB6')
two_ano_test_selected       <- stats_list$two_ano_test_selected
tukey_test_selected         <- stats_list$tukey_test_selected
Proportinal_plot(IntData, 'IntData', 'class_new', c(-4.5, 1), 20, 6, IntData_class_prop_filtered,two_ano_test_selected,tukey_test_selected, 'Prop_noB6', 'for_figure')
write.csv(IntData_class_prop_filtered, file = paste0(output, '/IntData_class_prop_filtered.csv'))

IntData_Th_supertype_prop <- Prorpotiona_analysis_create_df(IntData_Th_supertype, 'IntData_Th_supertype', 'marker') 
stats_list                <- Prorpotiona_analysis_stat(IntData_Th_supertype_prop$Prop_noB6, 'IntData_Th_supertype', my_comparisons, 'Prop_noB6')
two_ano_test_selected     <- stats_list$two_ano_test_selected
tukey_test_selected       <- stats_list$tukey_test_selected
Proportinal_plot(IntData_Th_supertype, 'IntData_Th_supertype', 'marker', c(-5.5, 2), 20, 6, IntData_Th_supertype_prop$Prop_noB6,two_ano_test_selected,tukey_test_selected, 'Prop_noB6', '')
df3 <- IntData_Th_supertype_prop$Prop_noB6
write.csv(df3, file = paste0(output, '/IntData_Th_supertype_prop.csv'))

IntData_Th_cluster_prop <- Prorpotiona_analysis_create_df(IntData_Th_cluster, 'IntData_Th_cluster', 'marker') 
stats_list                <- Prorpotiona_analysis_stat(IntData_Th_cluster_prop$Prop_noB6, 'IntData_Th_cluster', my_comparisons, 'Prop_noB6')
two_ano_test_selected     <- stats_list$two_ano_test_selected
tukey_test_selected       <- stats_list$tukey_test_selected
Proportinal_plot(IntData_Th_cluster, 'IntData_Th_cluster', 'marker', c(-6, -1), 20, 6, IntData_Th_cluster_prop$Prop_noB6,two_ano_test_selected,tukey_test_selected, 'Prop_noB6', '')
df4 <- IntData_Th_cluster_prop$Prop_noB6
write.csv(df4, file = paste0(output, '/IntData_Th_cluster_prop.csv'))

IntData_Th_supertype_prop_ctrl <- Prorpotiona_analysis_create_df(IntData_Th_supertype, 'IntData_Th_supertype', 'marker') 
stats_list                <- Prorpotiona_analysis_stat(IntData_Th_supertype_prop_ctrl$Prop_only_ctrls, 'IntData_Th_supertype', my_comparisons, 'Prop_only_ctrls')
two_ano_test_selected     <- stats_list$two_ano_test_selected
tukey_test_selected       <- stats_list$tukey_test_selected
Proportinal_plot(IntData_Th_supertype, 'IntData_Th_supertype', 'marker', c(-4.5, 0), 20, 6, IntData_Th_supertype_prop_ctrl$Prop_only_ctrls,two_ano_test_selected,tukey_test_selected, 'Prop_only_ctrls', '')
df5 <- IntData_Th_supertype_prop_ctrl$Prop_only_ctrls
write.csv(df5, file = paste0(output, '/IntData_Th_supertype_prop_ctrls.csv'))

IntData_Th_cluster_prop_ctrl <- Prorpotiona_analysis_create_df(IntData_Th_cluster, 'IntData_Th_cluster', 'marker') 
stats_list                <- Prorpotiona_analysis_stat(IntData_Th_cluster_prop_ctrl$Prop_only_ctrls, 'IntData_Th_cluster', my_comparisons, 'Prop_only_ctrls')
two_ano_test_selected     <- stats_list$two_ano_test_selected
tukey_test_selected       <- stats_list$tukey_test_selected
Proportinal_plot(IntData_Th_cluster, 'IntData_Th_cluster', 'marker',  c(-4.5, -1), 20, 6, IntData_Th_cluster_prop_ctrl$Prop_only_ctrls,two_ano_test_selected,tukey_test_selected, 'Prop_only_ctrls', '')
df6 <- IntData_Th_cluster_prop_ctrl$Prop_only_ctrls
write.csv(df6, file = paste0(output, '/IntData_Th_cluster_prop_ctrls.csv'))

# GO analysis using gprofiler: https://biit.cs.ut.ee/gprofiler/gost
pre_GO_analysis(IntData_Th_cluster, 'IntData_Th_cluster', 0.01)
GO_analysis_csv(IntData_Th_cluster, 'IntData_Th_cluster', 0.01)
GO_analysis_images(IntData_Th_cluster, 'IntData_Th_cluster', 0.01) 