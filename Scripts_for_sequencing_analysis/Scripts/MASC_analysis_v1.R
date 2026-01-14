MASC_updated <- function(dataset, cluster, contrast, random_effects, fixed_effects = NULL, verbose = FALSE, save_models = FALSE,
                  save_model_dir = NULL){
  if(is.factor(dataset[[contrast]]) == FALSE) {
    stop('Specified contrast term is not coded as a factor in dataset')
  }
  cluster    <- as.character(cluster)
  designmat  <- model.matrix(~cluster + 0, data.frame(cluster = cluster))
  dataset    <- cbind(designmat, dataset)
  cluster    <- as.character(cluster)
  designmat  <- model.matrix(~cluster + 0, data.frame(cluster = cluster))
  dataset    <- cbind(designmat, dataset)
  res        <- vector(mode = 'list', length = length(unique(cluster)))
  names(res) <- attributes(designmat)$dimnames[[2]]
  if(!is.null(fixed_effects) && !is.null(random_effects)){
     model_rhs <- paste0(c(paste0(fixed_effects, collapse = ' + '), paste0('(1|', random_effects, ')', collapse = ' + ')), collapse = ' + ')
    if(verbose == TRUE) {
      message(paste('Using null model:', 'cluster ~', model_rhs))
    }
  } else if(!is.null(fixed_effects) && is.null(random_effects)){
    model_rhs <- paste0(fixed_effects, collapse = ' + ')
    if (verbose == TRUE) {
      message(paste('Using null model:', 'cluster ~', model_rhs))
      stop('No random effects specified')
    }
  } else if(is.null(fixed_effects) && !is.null(random_effects)){
    model_rhs <- paste0('(1|', random_effects, ')', collapse = ' + ')
    if (verbose == TRUE) {
      message(paste('Using null model:', 'cluster ~', model_rhs))
    }
  } else{
    model_rhs <- '1'
    if(verbose == TRUE) {
      message(paste('Using null model:', 'cluster ~', model_rhs))
      stop('No random or fixed effects specified')
    }
  }
  cluster_models <- vector(mode = 'list', length = length(attributes(designmat)$dimnames[[2]]))
  names(cluster_models) <- attributes(designmat)$dimnames[[2]]
  print('Beginning Modelling')
  for(i in seq_along(attributes(designmat)$dimnames[[2]])) {
    test_cluster <- attributes(designmat)$dimnames[[2]][i]
    if(verbose == TRUE) {
      message(paste('Creating logistic mixed models for', test_cluster))
    }
    dataset[,i]   <- factor(dataset[,i])
    null_fm       <- as.formula(paste0(c(paste0(test_cluster, ' ~ 1 + '), model_rhs), collapse = ''))
    full_fm       <- as.formula(paste0(c(paste0(test_cluster, ' ~ ', contrast, ' + '), model_rhs), collapse = ''))
    null_model    <- lme4::glmer(formula = null_fm, data = dataset, family = binomial, nAGQ = 1, verbose = 0, control = glmerControl(optimizer = 'bobyqa'))
    full_model    <- lme4::glmer(formula = full_fm, data = dataset, family = binomial, nAGQ = 1, verbose = 0, control = glmerControl(optimizer = 'bobyqa'))
    model_lrt     <- anova(null_model, full_model)
    contrast_lvl2 <- paste0(contrast, levels(dataset[[contrast]])[2])
    contrast_ci   <- confint.merMod(full_model, method = 'Wald', parm = contrast_lvl2)
    cluster_models[[i]]$null_model <- null_model
    cluster_models[[i]]$full_model <- full_model
    cluster_models[[i]]$model_lrt  <- model_lrt
    cluster_models[[i]]$confint    <- contrast_ci
  }
  output <- data.frame(cluster = attributes(designmat)$dimnames[[2]], size = colSums(designmat))
  output$model.pvalue                                                    <- sapply(cluster_models, function(x) x$model_lrt[['Pr(>Chisq)']][2])
  output[[paste(contrast_lvl2, 'OR', sep = '.')]]                        <- sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))
  output[[paste(contrast_lvl2, 'OR', '95pct.ci.lower', sep = '.')]]      <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, '2.5 %']))
  output[[paste(contrast_lvl2, 'OR', '95pct.ci.upper', sep = '.')]]      <- sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, '97.5 %']))
  output[[paste(contrast_lvl2, 'log2_OR', sep = '.')]]                   <- log2(sapply(cluster_models, function(x) exp(fixef(x$full)[[contrast_lvl2]]))) # added: log scale
  output[[paste(contrast_lvl2, 'log2_OR', '95pct.ci.lower', sep = '.')]] <- log2(sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, '2.5 %']))) # added: log scale
  output[[paste(contrast_lvl2, 'log2_OR', '95pct.ci.upper', sep = '.')]] <- log2(sapply(cluster_models, function(x) exp(x$confint[contrast_lvl2, '97.5 %']))) # added: log scale
  output$cluster <- gsub('cluster','',as.character(output$cluster)) # added: remove string of 'cluster' in the output
  
  if(save_models == TRUE){
    saveModelObj(cluster_models, save_dir = save_model_dir)
    return(output)
  }
  else{
    return(output)
  }
}

Run_MASC <- function(data, groups_to_keep, output_file_suffix, annotation) {
  meta   <- data[data$groups %in% groups_to_keep, ] # Subset data
  meta$groups    <- factor(meta$groups)
  meta$Sex       <- factor(meta$Sex)
  meta$age       <- factor(meta$age)
  meta$sample_id <- factor(meta$sample_id)

  MASC_res <- MASC_updated(data = meta, cluster = meta[[annotation]], contrast = 'groups', random_effects = 'sample_id', verbose = TRUE, fixed_effects = c('Sex')) %>%
    mutate(fdr = p.adjust(model.pvalue, method = 'fdr')) %>%
    arrange(model.pvalue)
  write.csv(MASC_res, file = paste0(results, '/output/MASC_res_', output_file_suffix, '_', annotation, '.csv')) 
}

Bubble_plot_log <- function(data, name, y, xlim_log=4, angle=45, element_text_y=15) {
  data$model.pvalue <- as.numeric(data$model.pvalue)
  data[paste0('groups', y, '.log2_OR')] <- as.numeric(unlist(data[paste0('groups', y, '.log2_OR')]))
  
  ggplot(data, aes(x = data[,paste0('groups', y, '.log2_OR')], y = cluster, color = fdr)) +
    labs(title = paste0('Odds Ratio: ', name), x = paste0('groups', y, '.log2_OR'), y = 'Cluster',  color = 'FDR') + 
    geom_point(size=3, alpha = 1) + 
    geom_errorbar(aes(xmin = data[,paste0('groups', y, '.log2_OR.95pct.ci.lower')], xmax = data[,paste0('groups', y, '.log2_OR.95pct.ci.upper')]), width = 0.5) +
    geom_vline(xintercept=0, linetype='dashed', colour='black', linewidth=1) +
    scale_size_continuous(range = c(2, 10)) +  
    scale_color_gradient(low = 'blue', high = 'red',  limits=c(0,1),  breaks = c(0, 1)) + 
    theme_minimal() + 
    theme(plot.title = element_text(size = 15, face = 'bold'), 
          axis.text.x = element_text(size = 16, angle = 0, hjust = 0.5), axis.text.y = element_text(size = element_text_y),
          axis.title.x = element_text(size = 12), axis.title.y = element_blank(),
          legend.position = c(1, 1), legend.box.background = element_rect(fill = 'white', color = 'white'), legend.justification = c('right', 'top'), legend.margin = margin(7, 7, 7, 7),
          legend.key.size = unit(0.15,'in'), legend.key.height = unit(0.15,'in'), legend.key.width = unit(0.15,'in'), legend.text = element_text(size=15)) + 
    xlim(-xlim_log, xlim_log) +
    xlab(paste0(y, '.log2_OR'))
}

MASC_analysis <- function(obj, objname, annotation, xlim, xlim_log, angle = 45, width_bubble, height_bubble, element_text_y) {
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
  
  # Adding sex information; Sample ID to Sex mapping
  sex_mapping <- c(
    'RA004' = 'Male', 'RA012' = 'Male', 'RA017' = 'Female', 'RA018' = 'Female',
    'RA005' = 'Male', 'RA013' = 'Male', 'RA006' = 'Female', 'RA021' = 'Female',
    'RA016' = 'Male', 'RA020' = 'Male', 'RA015' = 'Female', 'RA019' = 'Female',
    'RA034' = 'Male', 'RA049' = 'Male', 'RA039' = 'Female', 'RA047' = 'Female',
    'RA040' = 'Male', 'RA043' = 'Male', 'RA041' = 'Female', 'RA048' = 'Female',
    'RA035' = 'Male', 'RA050' = 'Male', 'RA036' = 'Female', 'RA044' = 'Female',
    'RA026' = 'Male', 'RA030' = 'Male', 'RA027' = 'Female', 'RA031' = 'Female',
    'RA025' = 'Male', 'RA051' = 'Male', 'RA032' = 'Female', 'RA046' = 'Female',
    'RA022' = 'Male', 'RA028' = 'Male', 'RA023' = 'Female', 'RA029' = 'Female'
  )
  sample_ids <- seurat_meta_data$sample_id 
  sex <- map_chr(sample_ids, ~sex_mapping[.x]) 
  seurat_meta_data$Sex <- sex 

  # Replaced the annotation so the MASC output gives output in context of LM (something to do with alphabetical order); the word of 'Disease' makes MitoPark first
  seurat_meta_data$groups <- gsub('d_MP_lm_8wk', '1_LM_8wk', seurat_meta_data$groups)
  seurat_meta_data$groups <- gsub('a_MP_8wk', '4_Disease_MP_8wk', seurat_meta_data$groups)
  seurat_meta_data$groups <- gsub('e_MP_lm_16wk', '2_LM_16wk', seurat_meta_data$groups)
  seurat_meta_data$groups <- gsub('b_MP_16wk', '5_Disease_MP_16wk', seurat_meta_data$groups)
  seurat_meta_data$groups <- gsub('f_MP_lm_24wk', '3_LM_24wk', seurat_meta_data$groups)
  seurat_meta_data$groups <- gsub('c_MP_24wk', '6_Disease_MP_24wk', seurat_meta_data$groups)
  seurat_meta_data$groups <- gsub('g_WT_8wk', '7_WT_8wk', seurat_meta_data$groups)
  seurat_meta_data$groups <- gsub('h_WT_16wk', '8_WT_16wk', seurat_meta_data$groups)
  seurat_meta_data$groups <- gsub('i_WT_24wk', '9_WT_24wk', seurat_meta_data$groups)
  
  Run_MASC(seurat_meta_data, c('1_LM_8wk', '4_Disease_MP_8wk'), paste0(objname, '_1_LM_8wk_4_Disease_MP_8wk'), annotation)
  Run_MASC(seurat_meta_data, c('2_LM_16wk', '5_Disease_MP_16wk'), paste0(objname, '_2_LM_16wk_5_Disease_MP_16wk'), annotation)
  Run_MASC(seurat_meta_data, c('3_LM_24wk', '6_Disease_MP_24wk'), paste0(objname, '_3_LM_24wk_6_Disease_MP_24wk'), annotation)
  Run_MASC(seurat_meta_data, c('4_Disease_MP_8wk', '5_Disease_MP_16wk'), paste0(objname, '_4_Disease_MP_8wk_5_Disease_MP_16wk'), annotation)
  Run_MASC(seurat_meta_data, c('5_Disease_MP_16wk', '6_Disease_MP_24wk'), paste0(objname, '_5_Disease_MP_16wk_6_Disease_MP_24wk'), annotation)
  Run_MASC(seurat_meta_data, c('4_Disease_MP_8wk', '6_Disease_MP_24wk'), paste0(objname, '_4_Disease_MP_8wk_6_Disease_MP_24wk'), annotation)
  Run_MASC(seurat_meta_data, c('1_LM_8wk', '7_WT_8wk'), paste0(objname, '_1_LM_8wk_7_WT_8wk'), annotation)
  Run_MASC(seurat_meta_data, c('2_LM_16wk', '8_WT_16wk'), paste0(objname, '_2_LM_16wk_8_WT_16wk'), annotation)
  Run_MASC(seurat_meta_data, c('3_LM_24wk', '9_WT_24wk'), paste0(objname, '_3_LM_24wk_9_WT_24wk'), annotation)
  
  MASC_8wks_lm     <- read.csv(file.path(output, paste0('MASC_res_', objname, '_1_LM_8wk_4_Disease_MP_8wk', '_', annotation, '.csv')), row.names = 1)
  MASC_16wks_lm    <- read.csv(file.path(output, paste0('MASC_res_', objname, '_2_LM_16wk_5_Disease_MP_16wk', '_', annotation, '.csv')), row.names = 1)
  MASC_24wks_lm    <- read.csv(file.path(output, paste0('MASC_res_', objname, '_3_LM_24wk_6_Disease_MP_24wk', '_', annotation, '.csv')), row.names = 1)
  MASC_8wks_16wks  <- read.csv(file.path(output, paste0('MASC_res_', objname, '_4_Disease_MP_8wk_5_Disease_MP_16wk', '_', annotation, '.csv')), row.names = 1)
  MASC_16wks_24wks <- read.csv(file.path(output, paste0('MASC_res_', objname, '_5_Disease_MP_16wk_6_Disease_MP_24wk', '_', annotation, '.csv')), row.names = 1)
  MASC_8wks_24wks  <- read.csv(file.path(output, paste0('MASC_res_', objname, '_4_Disease_MP_8wk_6_Disease_MP_24wk', '_', annotation, '.csv')), row.names = 1)
  MASC_8wks_lm_WT  <- read.csv(file.path(output, paste0('MASC_res_', objname, '_1_LM_8wk_7_WT_8wk', '_', annotation, '.csv')), row.names = 1)
  MASC_16wks_lm_WT <- read.csv(file.path(output, paste0('MASC_res_', objname, '_2_LM_16wk_8_WT_16wk', '_', annotation, '.csv')), row.names = 1)
  MASC_24wks_lm_WT <- read.csv(file.path(output, paste0('MASC_res_', objname, '_3_LM_24wk_9_WT_24wk', '_', annotation, '.csv')), row.names = 1)
  
  png(file.path(output, paste0('MASC_', objname, '_', annotation, '_dotplots_log2_part1.png')), width = width_bubble, height = height_bubble/3, units = 'in', res = 300)
  p2 <- Bubble_plot_log(MASC_16wks_lm, 'MASC_16wks_lm','5_Disease_MP_16wk', xlim_log, element_text_y)
  p3 <- Bubble_plot_log(MASC_24wks_lm, 'MASC_24wks_lm','6_Disease_MP_24wk', xlim_log, element_text_y)
  p4 <- Bubble_plot_log(MASC_8wks_16wks, 'MASC_8wks_16wks','5_Disease_MP_16wk', xlim_log, element_text_y)
  print(p2+p3+p4)
  dev.off()
  
  png(file.path(output, paste0('MASC_', objname, '_', annotation, '_dotplots_log2_part2.png')), width = width_bubble, height = height_bubble/3, units = 'in', res = 300)
  p1 <- Bubble_plot_log(MASC_8wks_lm, 'MASC_8wks_lm', '4_Disease_MP_8wk', xlim_log, element_text_y)
  p2 <- Bubble_plot_log(MASC_16wks_lm, 'MASC_16wks_lm','5_Disease_MP_16wk', xlim_log, element_text_y)
  p3 <- Bubble_plot_log(MASC_24wks_lm, 'MASC_24wks_lm','6_Disease_MP_24wk', xlim_log, element_text_y)
  print(p1+p2+p3)
  dev.off()
  
  png(file.path(output, paste0('MASC_', objname, '_', annotation, '_dotplots_log2_part3.png')), width = width_bubble, height = height_bubble/3, units = 'in', res = 300)
  p7 <- Bubble_plot_log(MASC_8wks_lm_WT, 'MASC_8wks_lm_WT','7_WT_8wk', xlim_log, element_text_y)
  p8 <- Bubble_plot_log(MASC_16wks_lm_WT, 'MASC_16wks_lm_WT','8_WT_16wk', xlim_log, element_text_y)
  p9 <- Bubble_plot_log(MASC_24wks_lm_WT, 'MASC_24wks_lm_WT','9_WT_24wk', xlim_log, element_text_y)
  print(p7+p8+p9)
  dev.off()
}

