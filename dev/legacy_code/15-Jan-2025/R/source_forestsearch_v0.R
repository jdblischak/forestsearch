
# Remove virtual twins and plotting functions

source_fs_functions<-function(file_loc="/R/"){
files_toload<-c("forestsearch_functions_v0.R","subgroup_search_v0.R","subgroup_consistency_v0.R","forest_search_v0.R",
"grf_functions_v0.R","sim_aft_gbsg-mod4_v0.R","oc_analyses_m4_FS4_v0.R","oc_analyses_m4-FourCUTS_FS4_v0.R",
"forestsearch_bootstrap_functions_v0.R","forestsearch_bootstrap-parallel_v0.R","km_resampling_functions.R",
"get_FSdata_v0.R","bootstrap_summary_v0.R","subgroup_results_summary_v0.R","KMsurvdiff_functions_v0.R",
"forest_search_cross-validation_v0.R","forest_bootstrap_dofuture_v0.R","forest_bootstrap_doparallel_v0.R","cox_cubicspline_v0.R","sim_weibspline_functions_v0.R")

for(ff in 1:length(files_toload)){
source(c(paste0(file_loc,files_toload[ff])))
}
}
      
