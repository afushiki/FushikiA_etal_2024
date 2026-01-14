pre_GO_analysis <-  function(obj, objname, pval_th = 0.01) {
  Idents(obj) <- 'seurat_clusters'
  exportDE    <- function(clusterNumber, pval_thresh){
    fileName  <- paste0('cl',clusterNumber, '_th_subset_de_markers_padj_',pval_thresh, '_', objname, '.csv')
    FindMarkers(obj, ident.1 = clusterNumber) %>% # or only.pos = TRUE
      filter(p_val_adj < pval_thresh) %>% 
      rownames() %>%
      as.data.frame() %>%
      write.table(file=file.path(go_analysis_output, fileName), row.names=F, col.names='genes')
  }
  for (i in 0:(max(as.numeric(obj@meta.data$seurat_clusters))-1)){  # seurat_cluster starts from zero
    exportDE(i, pval_th)
  }
}

GO_analysis_csv <-  function(obj, objname, pval_th) {
  temp     <- list.files(go_analysis_output, pattern=paste0('_th_subset_de_markers_padj_',  pval_th,  '_', objname, '.csv'), full.names = T)
  filename <- gsub('.csv', '', basename(temp))
  
  temp_files <- list()
  for (i in 1:length(temp)){ temp_files[[i]] = read.csv(temp[i])  }
  
  go_result  <- list()
  go_result2 <- list()
  go_result_reordered <- list()
  for (i in 1:length(temp_files)){
    go_result[[i]] <- gost(query = temp_files[[i]]$genes, 
                           organism = 'mmusculus', ordered_query = FALSE, 
                           multi_query = FALSE, significant = TRUE, exclude_iea = FALSE, 
                           measure_underrepresentation = FALSE, evcodes = FALSE, 
                           user_threshold = 0.05, correction_method = 'g_SCS', 
                           domain_scope = 'annotated', custom_bg = NULL, 
                           numeric_ns = '', sources = NULL, as_short_link = FALSE, highlight = TRUE)
    go_result2[[i]]          <- as.data.frame(go_result[[i]]$result, stringsAsFactors = TRUE) 
    go_result_reordered[[i]] <- go_result2[[i]][, c('source',  'term_name', 'term_id',  'p_value')]
    go_result_reordered[[i]] <- go_result_reordered[[i]]  %>% mutate(negative_log10_of_adjusted_p_value = -log10(p_value))
    write.csv(go_result_reordered[[i]],  file.path(go_analysis_output, paste0(filename[i], '_gProfiler_results.csv')), row.names = FALSE)
  }
  
  temp_result <- list.files(go_analysis_output, pattern='_gProfiler_results.csv$', full.names = T)
  gProfiler_output_files <- list()

  for (i in 1:length(temp_result)){
    gProfiler_output_files[[i]] <- read.csv(temp_result[i]) %>% mutate(cluster = levels(factor(obj@meta.data$marker))[i])
  }

  gProfiler_output_results <- bind_rows(gProfiler_output_files)
  top_GO <- gProfiler_output_results %>%
    group_by(cluster,source) %>%
    slice_max(order_by = negative_log10_of_adjusted_p_value, n = 3) %>%
    select(1:3,cluster,negative_log10_of_adjusted_p_value)

  write.csv(gProfiler_output_results, file.path(go_analysis_output, paste0('gProfiler_all_results', '_', objname, '.csv')))  
  write.csv(top_GO, file.path(go_analysis_output, paste0('Top_gProfiler_final_results',  '_', objname, '.csv')))
  }

GO_analysis_images <-  function(obj, objname, pval_th) {
  top_GO <- read.csv(file.path(go_analysis_output, paste0('Top_gProfiler_final_results',  '_', objname, '.csv')))
  top_GO$negative_log10_of_adjusted_p_value <- format(round(top_GO$negative_log10_of_adjusted_p_value, 2), nsmall=2)
  colnames(top_GO) <- c('X', 'Source','Term name','Term ID', 'Cluster', 'log10(1/padj)') 
  top_GO <- top_GO[, !(names(top_GO) %in% 'X')]
  
  top_GO_split <- split(top_GO, top_GO$Cluster)
  gt_tbl_split <- list()
  
  for (i in 1:length(top_GO_split)){
    gt_tbl_split[[i]] <- gt(top_GO_split[[i]]) %>% tab_style(style = cell_text(weight = 'bold'), locations = cells_column_labels(columns=colnames(top_GO_split[[i]]))) %>% 
      tab_header(title = 'GO Table', subtitle = 'Based on dopamine subcluster') %>% opt_stylize(style = 1, color = 'gray') 
    gtsave(gt_tbl_split[[i]], file.path(go_analysis_output, paste0('GO_table_Th_subcluster_', names(top_GO_split)[i], '_', pval_th,  '_', objname, '.png')))
  }
}