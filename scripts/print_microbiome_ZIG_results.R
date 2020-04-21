##########################
## Statistical Analyses ##
##########################

for( a in 1:length(microbiome_statistical_analyses) ) {
  meg_fitZig(data_list=microbiome_analytic_data,
             data_names=microbiome_analytic_names,
             metadata=microbiome_metadata,
             zero_mod=model.matrix(~1 + log(libSize(microbiome))),
             data_mod=microbiome_statistical_analyses[[a]]$model_matrix,
             filter_min_threshold=0.15,
             contrast_list=microbiome_statistical_analyses[[a]]$contrasts,
             random_effect_var=microbiome_statistical_analyses[[a]]$random_effect,
             outdir=paste(stats_output_dir, 'Microbiome', microbiome_statistical_analyses[[a]]$name,
                          sep='/', collapse=''),
             analysis_name=microbiome_statistical_analyses[[a]]$name,
             analysis_subset=microbiome_statistical_analyses[[a]]$subsets,
             data_type='Microbiome',
             pval=0.99,
             top_hits=1000)
}
