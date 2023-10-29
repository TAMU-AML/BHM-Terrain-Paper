SOURCE_DIR = Sys.getenv("SOURCE_DIR")
LOG_DIR = Sys.getenv("LOG_DIR")
DATA_DIR = Sys.getenv("DATA_DIR")
ARTIFACTS_DIR = Sys.getenv("ARTIFACTS_DIR")
RESULTS_DIR = Sys.getenv("RESULTS_DIR")

library(latex2exp)
source(paste0(SOURCE_DIR, '/BHM-new-subroutines.R'))
terrain_variables = c("slope", "rix", "ridge")
file_name_dict = c(
  "slope_param_for_terrain_slope"="Figure8a.png",
  "location_param_for_terrain_slope"="Figure8b.png",
  "slope_param_for_terrain_rix"="Figure8c.png",
  "location_param_for_terrain_rix"="Figure8d.png",
  "slope_param_for_terrain_ridge"="Figure8e.png",
  "location_param_for_terrain_ridge"="Figure8f.png"
)
for (terrain in 1:length(terrain_variables)){
  betas = c()
  plot_labels = c("all")
  bhm_model = readRDS(paste0(ARTIFACTS_DIR, '/BHM-2dim-None-all.rds'))
  mc_index = get_mc_index(bhm_model) 
  betas_slope = bhm_model$samples$beta[terrain,1,mc_index]
  betas_location = bhm_model$samples$beta[terrain,2,mc_index]
  for (other in setdiff(c(1:length(terrain_variables)), terrain)){
    terrain_idx = sort(c(terrain, other))
    var_labels = paste(terrain_idx, collapse=',')
    plot_labels = c(plot_labels, paste(terrain_variables[terrain_idx], collapse=','))
    bhm_model = readRDS(paste0(ARTIFACTS_DIR, '/BHM-2dim-None-',var_labels,'.rds'))
    mc_index = get_mc_index(bhm_model) 
    betas_slope = cbind(betas_slope, bhm_model$samples$beta[which(terrain_idx==terrain),1,mc_index])
    betas_location = cbind(betas_location,bhm_model$samples$beta[which(terrain_idx==terrain),2,mc_index] )
  }
  bhm_model = readRDS(paste0(ARTIFACTS_DIR, '/BHM-2dim-None-',terrain,'.rds'))
  mc_index = get_mc_index(bhm_model) 
  betas_slope = cbind(betas_slope, bhm_model$samples$beta[1,1,mc_index])
  betas_location = cbind(betas_location,bhm_model$samples$beta[1,2,mc_index] )
  plot_labels = c(plot_labels, terrain_variables[terrain])
  
  png(paste0(RESULTS_DIR, '/', file_name_dict[paste0('slope_param_for_terrain_', terrain_variables[terrain])]), res=600, height=4, width=5.2, units='in')
  par(mai=c(0.85,1.4,0.75,0.85))
  boxplot(betas_slope, use.cols = TRUE, names=plot_labels, horizontal=TRUE, las=1, 
          xlab="Coefficient value", ylab="",
          main=bquote(paste("Terrain ", .(terrain_variables[terrain]),"'s coefficient for inflection slope, ", theta[1], sep="")), 
          ylim=c(-0.025,0.02))
  mtext("Terrain combination", side=2, line=5.5)
  dev.off()

  png(paste0(RESULTS_DIR, '/', file_name_dict[paste0('location_param_for_terrain_', terrain_variables[terrain])]), res=600, height=4, width=5.2, units='in')
  par(mai=c(0.85,1.4,0.75,0.85))
  boxplot(betas_location, use.cols = TRUE, names=plot_labels, horizontal=TRUE, 
          las=1, xlab="Coefficient value", 
          main=bquote(paste("Terrain ", .(terrain_variables[terrain]),"'s coefficient for inflection location, ", theta[2], sep="")),
          ylim=c(-0.2,0.3))
  mtext("Terrain combination", side=2, line=5.5)
  dev.off()
}
