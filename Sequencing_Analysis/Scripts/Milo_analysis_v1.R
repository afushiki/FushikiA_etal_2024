milo_analysis <-  function(obj, objname, reduction, groupby, threshold_frac, colours, pt.size = 2, label = TRUE, 
                           label_size, legend_ylim, legend_limit, legend_breaks, xaxis_size) {

  obj              <- Convert_Assay(seurat_object = obj, convert_to = 'V3')
  obj_sce          <- as.SingleCellExperiment(obj)
  obj_milo         <- Milo(obj_sce)
  obj_milo         <- buildGraph(obj_milo, k=30, d=30) 
  obj_milo         <- makeNhoods(obj_milo, prop=0.1, k=30, d=30, refined = TRUE) 
  obj_milo         <- countCells(obj_milo, sample = 'cell_ids', meta.data = obj@meta.data) 
  obj_design       <- data.frame(colData(obj_milo))[,c('cell_ids', 'groups')]
  obj_design       <- obj_design[, c('cell_ids', 'groups')]
  obj_design       <- distinct(obj_design); rownames(obj_design) <- obj_design$cell_ids
  obj_design       <- obj_design[colnames(nhoodCounts(obj_milo)), , drop=FALSE]
  obj_milo         <- calcNhoodDistance(obj_milo, d=30)
  da_results       <- testNhoods(obj_milo, design = ~ groups, design.df = obj_design)
  obj_milo         <- buildNhoodGraph(obj_milo)
  da_results       <- annotateNhoods(obj_milo, da_results, coldata_col = 'ident')
  da_results %>% arrange(FDR) %>% write_csv(file.path(output, paste0('da_results_', objname, '.csv'))) 
  da_results       <- da_results %>% arrange(desc(da_results$ident)) 
  da_results       <- da_results %>% mutate(across(-ident, as.numeric))
  da_results_copy  <- da_results # for subset
  da_results$ident <- ifelse(da_results$ident_fraction < threshold_frac, 'Mixed', da_results$ident)
  da_results$ident <- as.factor(da_results$ident)
  
  png(file.path(output, paste0('milo_', objname, '.png')), width = 18, height = 6, units = 'in', res = 300)
  image1 <- DimPlot(obj, reduction = reduction, group.by=groupby, split.by = 'groups', pt.size = pt.size, label=label, label.size = label_size, cols=colours, repel = TRUE) + 
            ggtitle(NULL) & NoLegend() & NoAxes()
  image2 <- plotNhoodGraphDA(obj_milo, da_results, alpha=0.1, size_range = c(2,6), node_stroke = 1) + 
            scale_fill_gradientn(name='logFC', colours = rev(brewer.pal(11, 'RdBu')),limits = legend_limit, breaks = legend_breaks) + 
            theme(legend.key.size = unit(0.75, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=20)) +
            guides(size=guide_legend(title='Nhood size',order = 1), edge_width=guide_legend(title='overlap size',order = 2)) 
  print(image1 + image2 + plot_layout(ncol=2, widths = c(2, 1)))
  dev.off()

  png(file.path(output, paste0('milo_beeswarm_', objname, '.png')), width = 12, height = 12, units = 'in', res = 300)
  print(plotDAbeeswarm(da_results, group.by = 'ident', alpha = 0.1) + 
        ggtitle(objname) + ylim(legend_ylim) + geom_point(aes(size = 6)) + 
        theme(axis.title.y=element_blank(), axis.text=element_text(size=xaxis_size), panel.grid.minor.y = element_blank()) +
        scale_color_gradient2(name='logFC', low = rev(brewer.pal(11, 'RdBu'))[1], 
                              mid= rev(brewer.pal(11, 'RdBu'))[6], high= rev(brewer.pal(11, 'RdBu'))[11], 
                              na.value = 'grey80', limits =legend_limit, breaks = legend_breaks)+ NoLegend())
  dev.off()
  
  if (groupby == 'cluster_name') {
    # subsetting for figures in 0882 supertype in ABC atlas
    da_results_0882  <- da_results_copy %>% filter(str_detect(ident, '3853_|3854_|3855_|3856_|3857_|3858_|3859_')) 
    da_results_0882$ident <- ifelse(da_results_0882$ident_fraction < threshold_frac, 'Mixed', da_results_0882$ident)
    da_results_0882$ident <- as.factor(da_results_0882$ident)
    
    print ("Creating subset figures...")
    png(file.path(output, paste0('milo_beeswarm_', objname, '_in_0882.png')), width = 8, height = 8, units = 'in', res = 300)
    print(plotDAbeeswarm(da_results_0882, group.by = 'ident', alpha = 0.1) + 
          ggtitle(objname) + ylim(legend_ylim) + geom_point(aes(size = 10)) +
          theme(axis.title.y=element_blank(), axis.text=element_text(size=xaxis_size+2), panel.grid.minor.y = element_blank()) +
          scale_color_gradient2(name='logFC', low = rev(brewer.pal(11, "RdBu"))[1],
                                mid= rev(brewer.pal(11, "RdBu"))[6], high= rev(brewer.pal(11, "RdBu"))[11], 
                                na.value = "grey80", limits =legend_limit, breaks = legend_breaks)+ NoLegend())
    dev.off()
    
    # since the number of clusters is large, use fewer colors and emphasize clusters with a fold change below −2
    groups_to_highlight <- unique(da_results_copy$ident[da_results_copy$logFC < -2])
    color_map           <- setNames(color_pal2, sort(unique(obj@meta.data$cluster_name))) 
    assigned_colors     <- ifelse(sort(unique(obj@meta.data$cluster_name)) %in% groups_to_highlight, color_map, '#cccccc') 
    non_grey            <- assigned_colors[assigned_colors != '#cccccc']
    grey                <- assigned_colors[assigned_colors == '#cccccc']
    assigned_colors_for_legend <- c(non_grey, grey)
    
    obj$custom_labels   <- Idents(obj)
    new_label           <- ifelse(levels(obj$custom_labels) %in% groups_to_highlight, levels(obj$custom_labels), "others")
    levels(obj$custom_labels) <- new_label
    assigned_colors_update2   <- c(non_grey, unique(grey)) # grey for others in label
    
    png(file.path(output, paste0('milo_', objname, '_label_limited_legend.png')), width = 21, height = 6, units = 'in', res = 300)
    image1 <- DimPlot(obj, reduction = reduction, group.by='custom_labels', split.by = 'groups', pt.size = pt.size, label=FALSE,
                      label.size = label_size +2 , cols=assigned_colors_update2, repel = TRUE, alpha = 0.6) +
              scale_color_manual(breaks = sort(levels(obj$custom_labels)), values = assigned_colors_update2) +
              theme(legend.text=element_text(size=18), legend.key.size = unit(1, 'line'), legend.position="right") + guides(color=guide_legend(ncol =1, override.aes = list(size = 6))) +
              ggtitle(NULL) & NoAxes()
    image2 <- plotNhoodGraphDA(obj_milo, da_results, alpha=0.1, size_range = c(2,6), node_stroke = 1) + 
              scale_fill_gradientn(name='logFC', colours = rev(brewer.pal(11, 'RdBu')),limits = legend_limit, breaks = legend_breaks) + 
              theme(legend.key.size = unit(0.75, 'cm'), legend.title = element_text(size=20), legend.text = element_text(size=20)) +
              guides(size=guide_legend(title='Nhood size',order = 1), edge_width=guide_legend(title='overlap size',order = 2)) 
    print(image1 + image2 + plot_layout(ncol=2, widths = c(2, 1)))
    dev.off()  
    
  } else {
    print ("No subsetting, done!")
  }
  
}