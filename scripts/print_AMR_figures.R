## Print out some exploratory figures 


#############################################
## Exploratory Analyses: Alpha Rarefaction ##
#############################################
for( v in 1:length(AMR_exploratory_analyses) ) {
  # AMR
  meg_alpha_rarefaction(data_list=AMR_raw_analytic_data,
                        data_names=AMR_raw_analytic_names,
                        metadata=metadata,
                        sample_var=sample_column_id,
                        group_var=AMR_exploratory_analyses[[v]]$exploratory_var,
                        analysis_subset=AMR_exploratory_analyses[[v]]$subsets,
                        outdir=paste(graph_output_dir, 'AMR', AMR_exploratory_analyses[[v]]$name,
                                     sep='/', collapse=''),
                        data_type='AMR')
}

######################################
## Exploratory Analyses: Ordination ##
######################################
for( v in 1:length(AMR_exploratory_analyses) ) {
  # AMR NMDS
  meg_ordination(data_list = AMR_analytic_data,
                 data_names = AMR_analytic_names,
                 metadata = metadata,
                 sample_var = sample_column_id,
                 hull_var = AMR_exploratory_analyses[[v]]$exploratory_var,
                 analysis_subset=AMR_exploratory_analyses[[v]]$subsets,
                 outdir = paste(graph_output_dir, 'AMR', AMR_exploratory_analyses[[v]]$name,
                                sep='/', collapse=''),
                 data_type = 'AMR',
                 method = 'NMDS')
  
  # AMR PCA
  meg_ordination(data_list = AMR_analytic_data,
                 data_names = AMR_analytic_names,
                 metadata = metadata,
                 sample_var = sample_column_id,
                 hull_var = AMR_exploratory_analyses[[v]]$exploratory_var,
                 analysis_subset=AMR_exploratory_analyses[[v]]$subsets,
                 outdir = paste(graph_output_dir, 'AMR', AMR_exploratory_analyses[[v]]$name,
                                sep='/', collapse=''),
                 data_type = 'AMR',
                 method = 'PCA')
}

####################################
## Exploratory Analyses: Heatmaps ##
####################################

# AMR Heatmaps for each level
for( v in 1:length(AMR_exploratory_analyses) ) {
  for( l in 1:length(AMR_analytic_names) ) {
    meg_heatmap(melted_data=amr_melted_analytic,
                metadata=metadata,
                sample_var=sample_column_id,
                group_var=AMR_exploratory_analyses[[v]]$exploratory_var,
                level_var=AMR_analytic_names[l],
                analysis_subset=AMR_exploratory_analyses[[v]]$subsets,
                outdir=paste(graph_output_dir, 'AMR', AMR_exploratory_analyses[[v]]$name,
                             sep='/', collapse=''),
                data_type='AMR')
  }
}



####################################
## Exploratory Analyses: Barplots ##
####################################

# AMR
for( v in 1:length(AMR_exploratory_analyses) ) {
  for( l in 1:length(AMR_analytic_names) ) {
    suppressWarnings(
      meg_barplot(melted_data=amr_melted_analytic,
                  metadata=metadata,
                  sample_var=sample_column_id,
                  group_var=AMR_exploratory_analyses[[v]]$exploratory_var,
                  level_var=AMR_analytic_names[l],
                  analysis_subset=AMR_exploratory_analyses[[v]]$subsets,
                  outdir=paste(graph_output_dir, 'AMR', AMR_exploratory_analyses[[v]]$name,
                               sep='/', collapse=''),
                  data_type='AMR')
    )
  }
}

##########################
## Statistical Analyses ##
##########################
for( a in 1:length(AMR_statistical_analyses) ) {
  meg_fitZig(data_list=AMR_analytic_data,
             data_names=AMR_analytic_names,
             metadata=metadata,
             zero_mod=model.matrix(~1 + log(libSize(amr))),
             data_mod=AMR_statistical_analyses[[a]]$model_matrix,
             filter_min_threshold=0.15,
             contrast_list=AMR_statistical_analyses[[a]]$contrasts,
             random_effect_var=AMR_statistical_analyses[[a]]$random_effect,
             outdir=paste(stats_output_dir, 'AMR', AMR_statistical_analyses[[a]]$name,
                          sep='/', collapse=''),
             analysis_name=AMR_statistical_analyses[[a]]$name,
             analysis_subset=AMR_statistical_analyses[[a]]$subsets,
             data_type='AMR',
             pval=0.99,
             top_hits=1000)
}


########################
## Output of matrices ##
########################
write.csv(make_sparse(amr_class, 'class', c('class')), 'amr_matrices/sparse_normalized/AMR_Class_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_class, 'amr_matrices/normalized/AMR_Class_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_class_raw, 'amr_matrices/raw/AMR_Class_Raw.csv', sep=',', row.names = F, col.names = T)


write.csv(make_sparse(amr_mech, 'mechanism', c('mechanism')), 'amr_matrices/sparse_normalized/AMR_Mechanism_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_mech, 'amr_matrices/normalized/AMR_Mechanism_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_mech_raw, 'amr_matrices/raw/AMR_Mechanism_Raw.csv', sep=',', row.names = F, col.names = T)

write.csv(make_sparse(amr_group, 'group', c('group')), 'amr_matrices/sparse_normalized/AMR_Group_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_group, 'amr_matrices/normalized/AMR_Group_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_mech_raw, 'amr_matrices/raw/AMR_Group_Raw.csv', sep=',', row.names = F, col.names = T)

write.csv(make_sparse(amr_norm, 'header', c('header', 'class', 'mechanism', 'group')),
          'amr_matrices/sparse_normalized/AMR_Gene_Sparse_Normalized.csv',
          row.names=T)
write.table(amr_norm, 'amr_matrices/normalized/AMR_Gene_Normalized.csv', sep=',', row.names = F, col.names = T)
write.table(amr_raw, 'amr_matrices/raw/AMR_Gene_Raw.csv', sep=',', row.names = F, col.names = T)

