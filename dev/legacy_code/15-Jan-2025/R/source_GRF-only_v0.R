
source_fs_functions<-function(file_loc="/R/"){
files_toload<-c("grf_functions_v0.R","km_resampling_functions.R","subgroup_results_summary_v0.R","plotting_functions_v0.R")

for(ff in 1:length(files_toload)){
source(c(paste0(file_loc,files_toload[ff])))
}
}
      