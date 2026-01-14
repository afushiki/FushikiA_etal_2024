Proportinal_plot <- function(obj, objname, annotation, legend_ylim_proportion, width_proportion, heigh_proportion, 
                                  prop_matrix,two_ano_test_selected,tukey_test_selected, animal_groups, note) {
  
  tukey_test_selected$y.position <- tukey_test_selected$y.position - 1.3
  
  if (animal_groups == 'Prop_noB6') {
    png(file.path(output, paste0('Proportion_plots_', objname, '_', annotation, '_w_statistics_', animal_groups, '_', note, '.png')), width = width_proportion, height = heigh_proportion-1, units = 'in', res = 300)
    p1 <- ggboxplot(prop_matrix, x='groups', y='log(prop)',  facet.by='cluster', fill='groups') +
      facet_wrap(~cluster, ncol=length(levels(prop_matrix$cluster))) + 
      stat_compare_means(aes(label = paste0("ANOVA p=", after_stat(two_ano_test_selected$p))), label.x.npc=0.20, size=5, label.y=legend_ylim_proportion[1]) +
      stat_pvalue_manual(tukey_test_selected, hide.ns = TRUE, label='p.adj.signif', label.size = 8, bracket.size = 0.8, tip.length = 0.05, step.increase = 0.1, step.group.by = "cluster") +
      scale_fill_manual(values = c('#FF8888', '#FB5A5A', '#F43030', '#86CEFA', '#73B9EE', '#5494DA')) +
      theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 20), axis.title.y=element_text(size=25), axis.title.x=element_blank(), 
            legend.key.size = unit(0.5,'in'), legend.key.height = unit(0.5,'in'), legend.key.width = unit(0.5,'in'), legend.text = element_text(size=15), legend.position='right', legend.title = element_blank(),
            strip.text = element_text(face='bold', size=15)) +
      guides(fill = guide_legend(nrow = 9)) +
      ylab('Proportion (log10)') +
      coord_cartesian(ylim = legend_ylim_proportion)
    print(p1)
    dev.off()
  } else {
    png(file.path(output, paste0('Proportion_plots_', objname, '_', annotation, '_w_statistics_', animal_groups, '_', note, '.png')), width = width_proportion, height = heigh_proportion-1, units = 'in', res = 300)
    p1 <- ggboxplot(prop_matrix, x='groups', y='log(prop)',  facet.by='cluster', fill='groups') +
      facet_wrap(~cluster, ncol=length(levels(prop_matrix$cluster))) + 
      stat_compare_means(aes(label = paste0("ANOVA p=", after_stat(two_ano_test_selected$p))), label.x.npc=0.20, size=5, label.y=legend_ylim_proportion[1]) +
      stat_pvalue_manual(tukey_test_selected, hide.ns = TRUE, label='p.adj.signif', label.size = 8, bracket.size = 0.8, tip.length = 0.05, step.increase = 0.1, step.group.by = "cluster") +
      scale_fill_manual(values = c('#86CEFA', '#73B9EE', '#5494DA', '#CDCFCC', '#B7B8B3', '#9A9C99')) +
      theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 20), axis.title.y=element_text(size=25), axis.title.x=element_blank(), 
            legend.key.size = unit(0.5,'in'), legend.key.height = unit(0.5,'in'), legend.key.width = unit(0.5,'in'), legend.text = element_text(size=15), legend.position='right', legend.title = element_blank(),
            strip.text = element_text(face='bold', size=15)) +
      guides(fill = guide_legend(nrow = 9)) +
      ylab('Proportion (log10)') +
      coord_cartesian(ylim = legend_ylim_proportion)      
    print(p1)
    dev.off()
  }
  
  if (objname == 'IntData_Th_cluster') {
    
    tukey_test_selected$y.position <- tukey_test_selected$y.position - 1.3
    
     if (animal_groups == 'Prop_noB6') {
       png(file.path(output, paste0('Proportion_plots_', objname, '_', annotation, '_w_statistics_', animal_groups, '_', note, '.png')), width = width_proportion, height = heigh_proportion-1, units = 'in', res = 300)
       p1 <- ggboxplot(prop_matrix, x='groups', y='log(prop)',  facet.by='cluster', fill='groups') +
         facet_wrap(~cluster, ncol=length(levels(prop_matrix$cluster))) + 
         stat_compare_means(aes(label = paste0("ANOVA p=", after_stat(two_ano_test_selected$p))), label.x.npc=0.25, size=5, label.y=legend_ylim_proportion[1]) +
         stat_pvalue_manual(tukey_test_selected, hide.ns = TRUE,  label='p.adj.signif', label.size = 8, bracket.size = 0.8, tip.length = 0.05, step.increase = 0.1, step.group.by = "cluster") +
         scale_fill_manual(values = c('#FF8888', '#FB5A5A', '#F43030', '#86CEFA', '#73B9EE', '#5494DA')) +
         theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 20), axis.title.y=element_text(size=25), axis.title.x=element_blank(), 
               legend.key.size = unit(0.5,'in'), legend.key.height = unit(0.5,'in'), legend.key.width = unit(0.5,'in'), legend.text = element_text(size=15), legend.position='right', legend.title = element_blank(),
               strip.text = element_text(face='bold', size=15)) +
         guides(fill = guide_legend(nrow = 9)) +
         ylab('Proportion (log10)') +
         coord_cartesian(ylim = legend_ylim_proportion)
       print(p1)
       dev.off()
     } else {
       png(file.path(output, paste0('Proportion_plots_', objname, '_', annotation, '_w_statistics_', animal_groups, '_', note, '.png')), width = width_proportion, height = heigh_proportion-1, units = 'in', res = 300)
       p1 <- ggboxplot(prop_matrix, x='groups', y='log(prop)',  facet.by='cluster', fill='groups') +
         facet_wrap(~cluster, ncol=length(levels(prop_matrix$cluster))) + 
         stat_compare_means(aes(label = paste0("ANOVA p=", after_stat(two_ano_test_selected$p))), label.x.npc=0.25, size=5, label.y=legend_ylim_proportion[1]) +
         stat_pvalue_manual(tukey_test_selected, hide.ns = TRUE,  label='p.adj.signif', label.size = 8, bracket.size = 0.8, tip.length = 0.05, step.increase = 0.1, step.group.by = "cluster") +
         scale_fill_manual(values = c('#86CEFA', '#73B9EE', '#5494DA', '#CDCFCC', '#B7B8B3', '#9A9C99')) +
         theme(axis.text.x = element_blank(), axis.text.y = element_text(size = 20), axis.title.y=element_text(size=25), axis.title.x=element_blank(), 
               legend.key.size = unit(0.5,'in'), legend.key.height = unit(0.5,'in'), legend.key.width = unit(0.5,'in'), legend.text = element_text(size=15), legend.position='right', legend.title = element_blank(),
               strip.text = element_text(face='bold', size=15)) +
         guides(fill = guide_legend(nrow = 9)) +
         ylab('Proportion (log10)') +
         coord_cartesian(ylim = legend_ylim_proportion)
       print(p1)
       dev.off()
     }
    print ("Subsetted, done!")
  } else {
    print ("No subsetting, done!")
  }
}

