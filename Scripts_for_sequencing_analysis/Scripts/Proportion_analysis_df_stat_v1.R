Prorpotiona_analysis_create_df <- function(obj, objname, annotation) {
  obj$age <- obj$groups
  obj$age <- gsub('a_MP_8wk', '8wk', obj$age)
  obj$age <- gsub('b_MP_16wk', '16wk', obj$age)
  obj$age <- gsub('c_MP_24wk', '24wk', obj$age)
  obj$age <- gsub('d_MP_lm_8wk', '8wk', obj$age)
  obj$age <- gsub('e_MP_lm_16wk', '16wk', obj$age)
  obj$age <- gsub('f_MP_lm_24wk', '24wk', obj$age)
  obj$age <- gsub('g_WT_8wk', '8wk', obj$age)
  obj$age <- gsub('h_WT_16wk', '16wk', obj$age)
  obj$age <- gsub('i_WT_24wk', '24wk', obj$age)
  seurat_meta_data <- obj@meta.data
  
  # Extracting sample IDs from column names
  sample_ids <- rownames(seurat_meta_data) %>% str_split('_') %>% map_chr(~ifelse(length(.) > 1, .[[1]], NA))
  # Creating a new column in the metadata data frame for sample IDs
  seurat_meta_data <- seurat_meta_data %>% mutate(sample_id = sample_ids)
  # Make a table with frequency data by region and sample
  boxplot_table <- table(seurat_meta_data[['sample_id']], seurat_meta_data[[annotation]], seurat_meta_data[['groups']])
  # Convert the frequency data table into a matrix to make a boxplot
  boxplot_df<-as.data.frame(boxplot_table); colnames(boxplot_df) <- c('sample', 'cluster', 'groups', 'frequency')
  
  SNboxplotmatrix_prop <- boxplot_df %>% group_by(sample, groups) %>%　mutate(prop = frequency / sum(frequency)) %>%　ungroup()
  
  # Rename for figures
  SNboxplotmatrix_prop$groups  <- gsub('a_MP_8wk', 'MitoPark 8wks', SNboxplotmatrix_prop$groups)
  SNboxplotmatrix_prop$groups  <- gsub('b_MP_16wk', 'MitoPark 16wks', SNboxplotmatrix_prop$groups)
  SNboxplotmatrix_prop$groups  <- gsub('c_MP_24wk', 'MitoPark 24wks', SNboxplotmatrix_prop$groups)
  SNboxplotmatrix_prop$groups  <- gsub('d_MP_lm_8wk', 'Littermate 8wks', SNboxplotmatrix_prop$groups)
  SNboxplotmatrix_prop$groups  <- gsub('e_MP_lm_16wk', 'Littermate 16wks', SNboxplotmatrix_prop$groups)
  SNboxplotmatrix_prop$groups  <- gsub('f_MP_lm_24wk', 'Littermate 24wks', SNboxplotmatrix_prop$groups)
  SNboxplotmatrix_prop$groups  <- gsub('g_WT_8wk', 'B6 8wks', SNboxplotmatrix_prop$groups)
  SNboxplotmatrix_prop$groups  <- gsub('h_WT_16wk', 'B6 16wks', SNboxplotmatrix_prop$groups)
  SNboxplotmatrix_prop$groups  <- gsub('i_WT_24wk', 'B6 24wks', SNboxplotmatrix_prop$groups)
  SNboxplotmatrix_prop$groups  <- factor(SNboxplotmatrix_prop$groups, levels = c('MitoPark 8wks', 'MitoPark 16wks', 'MitoPark 24wks', 
                                                                                 'Littermate 8wks', 'Littermate 16wks' ,'Littermate 24wks', 
                                                                                 'B6 8wks', 'B6 16wks', 'B6 24wks'))
  SNboxplotmatrix_prop$cluster <- factor(SNboxplotmatrix_prop$cluster, levels = sort(levels(SNboxplotmatrix_prop$cluster)))
  SNboxplotmatrix_prop         <- filter(SNboxplotmatrix_prop, prop>0)
  SNboxplotmatrix_prop_rebuild <- SNboxplotmatrix_prop %>% mutate(genotype = case_when(grepl('MitoPark', groups) ~ 'MitoPark', grepl('Littermate', groups) ~ 'Control1', grepl('B6', groups) ~ 'Control2'))
  SNboxplotmatrix_prop_rebuild <- SNboxplotmatrix_prop_rebuild %>% mutate(age = case_when(grepl('8wks', groups) ~ '8wks', grepl('16wks', groups) ~ '16wks',  grepl('24wks', groups) ~ '24wks'))
  
  # Remove B6 animals for main figure
  SNboxplotmatrix_prop_noB6         <- SNboxplotmatrix_prop[!grepl('B6', SNboxplotmatrix_prop$groups),] 
  SNboxplotmatrix_prop_noB6_rebuild <- SNboxplotmatrix_prop_noB6 %>% mutate(genotype = case_when(grepl('MitoPark', groups) ~ 'MitoPark', grepl('Littermate', groups) ~ 'Control'))
  SNboxplotmatrix_prop_noB6_rebuild <- SNboxplotmatrix_prop_noB6_rebuild %>% mutate(age = case_when(grepl('8wks', groups) ~ '8wks', grepl('16wks', groups) ~ '16wks',  grepl('24wks', groups) ~ '24wks'))
  
  # Remove MitoPark animals for supplemental figure
  SNboxplotmatrix_prop_only_ctrls         <- SNboxplotmatrix_prop[!grepl('MitoPark', SNboxplotmatrix_prop$groups),] 
  SNboxplotmatrix_prop_only_ctrls_rebuild <- SNboxplotmatrix_prop_only_ctrls %>% mutate(genotype = case_when(grepl('Littermate', groups) ~ 'Control1', grepl('B6', groups) ~ 'Control2'))
  SNboxplotmatrix_prop_only_ctrls_rebuild <- SNboxplotmatrix_prop_only_ctrls_rebuild %>% mutate(age = case_when(grepl('8wks', groups) ~ '8wks', grepl('16wks', groups) ~ '16wks',  grepl('24wks', groups) ~ '24wks'))
  
  if (objname == 'IntData_Th_cluster') {
    # subsetting for figures in 0882 sypertype in ABC atlas
    SNboxplotmatrix_prop_copy  <- SNboxplotmatrix_prop
    SNboxplotmatrix_prop_0882  <- SNboxplotmatrix_prop_copy %>% filter(str_detect(cluster, '3853_|3854_|3855_|3856_|3857_|3858_|3859_')) 
    SNboxplotmatrix_prop_0882$cluster <- as.factor(SNboxplotmatrix_prop_0882$cluster)
    SNboxplotmatrix_prop_0882         <- filter(SNboxplotmatrix_prop_0882, prop>0)
    SNboxplotmatrix_prop_0882_rebuild <- SNboxplotmatrix_prop_0882 %>% mutate(genotype = case_when(grepl('MitoPark', groups) ~ 'MitoPark', grepl('Littermate', groups) ~ 'Control1', grepl('B6', groups) ~ 'Control2'))
    SNboxplotmatrix_prop_0882_rebuild <- SNboxplotmatrix_prop_0882_rebuild %>% mutate(age = case_when(grepl('8wks', groups) ~ '8wks', grepl('16wks', groups) ~ '16wks',  grepl('24wks', groups) ~ '24wks'))
    
    # Remove B6 animals for main figure
    SNboxplotmatrix_prop_0882_noB6         <- SNboxplotmatrix_prop_0882[!grepl("B6", SNboxplotmatrix_prop_0882$groups),]
    SNboxplotmatrix_prop_0882_noB6_rebuild <- SNboxplotmatrix_prop_0882_noB6 %>% mutate(genotype = case_when(grepl('MitoPark', groups) ~ 'MitoPark', grepl('Littermate', groups) ~ 'Control'))
    SNboxplotmatrix_prop_0882_noB6_rebuild <- SNboxplotmatrix_prop_0882_noB6_rebuild %>% mutate(age = case_when(grepl('8wks', groups) ~ '8wks', grepl('16wks', groups) ~ '16wks',  grepl('24wks', groups) ~ '24wks'))
    
    # Remove MitoPark animals for supplemental figure
    SNboxplotmatrix_prop_0882_only_ctrls         <- SNboxplotmatrix_prop_0882[!grepl('MitoPark', SNboxplotmatrix_prop_0882$groups),] 
    SNboxplotmatrix_prop_0882_only_ctrls_rebuild <- SNboxplotmatrix_prop_0882_only_ctrls %>% mutate(genotype = case_when(grepl('Littermate', groups) ~ 'Control1', grepl('B6', groups) ~ 'Control2'))
    SNboxplotmatrix_prop_0882_only_ctrls_rebuild <- SNboxplotmatrix_prop_0882_only_ctrls_rebuild %>% mutate(age = case_when(grepl('8wks', groups) ~ '8wks', grepl('16wks', groups) ~ '16wks',  grepl('24wks', groups) ~ '24wks'))
    
    print ("Subsetted, done!")
    return_list <- list()
    return_list <- list('Prop_noB6' = SNboxplotmatrix_prop_0882_noB6_rebuild, 'Prop_only_ctrls' = SNboxplotmatrix_prop_0882_only_ctrls_rebuild)
    return(return_list)
    
  } else {
    
    print ("No subsetting, done!")
    return_list <- list()
    return_list <- list('Prop_noB6' = SNboxplotmatrix_prop_noB6_rebuild, 'Prop_only_ctrls' = SNboxplotmatrix_prop_only_ctrls_rebuild)
    return(return_list)
    
  }
}

Prorpotiona_analysis_stat <- function(df, annotation, my_comparisons, animal_groups) {
  # Anova test
  two_ano_test_all      <- df %>% group_by(cluster) %>% anova_test(prop ~ genotype*age) %>% add_significance() # two-way ANOVA
  two_ano_test_selected <- subset(two_ano_test_all, Effect %in% 'genotype:age')
  ano_significant       <- subset(two_ano_test_selected, p < 0.05)
  # Tukey test
  tukey_test_all        <- df %>% group_by(cluster) %>% tukey_hsd(prop ~ groups) %>% add_xy_position()
  df_for_select         <- as.data.frame(my_comparisons, col.names = c('group1', 'group2'))
  tukey_test_selected   <- right_join(tukey_test_all, df_for_select, by = c('group1', 'group2')) # pick up only my comparisons
  tukey_test_selected   <- right_join(tukey_test_selected, ano_significant, by = 'cluster') # pick up only ANOVA significant

  two_ano_test_all_for_csv <- data.frame(lapply(two_ano_test_all, as.character), stringsAsFactors=FALSE)
  write.csv(two_ano_test_all_for_csv, file = paste0(output, '/two_ano_test_all_proportional_plots', annotation, '_', animal_groups, '.csv'))
  
  tukey_test_all_for_csv   <- data.frame(lapply(tukey_test_all, as.character), stringsAsFactors=FALSE)
  write.csv(tukey_test_all_for_csv, file = paste0(output, '/tukey_test_all_proportional_plots', annotation, '_', animal_groups, '.csv'))
  
  return_list <- list()
  return_list <- list('two_ano_test_selected' = two_ano_test_selected, 'tukey_test_selected' = tukey_test_selected)
  return(return_list)
}