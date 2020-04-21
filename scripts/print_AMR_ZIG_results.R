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
