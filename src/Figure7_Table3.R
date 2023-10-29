SOURCE_DIR = Sys.getenv("SOURCE_DIR")
LOG_DIR = Sys.getenv("LOG_DIR")
DATA_DIR = Sys.getenv("DATA_DIR")
ARTIFACTS_DIR = Sys.getenv("ARTIFACTS_DIR")
RESULTS_DIR = Sys.getenv("RESULTS_DIR")

source(paste0(SOURCE_DIR, '/BHM-new-subroutines.R'))
terrain_variables = c("slope", "rix", "ridge")
bhm_model = readRDS(paste0(ARTIFACTS_DIR, '/BHM-2dim-None-all.rds'))
mc_index = get_mc_index(bhm_model) 
beta = bhm_model$samples$beta
write.table(round(t(apply(beta[,,mc_index], c(1,2), mean)),6), file=paste0(RESULTS_DIR, "/Table3.txt"),
            sep='\t', row.names = c("beta_1", "beta_2"), col.names = terrain_variables)

png(paste0(RESULTS_DIR, '/Figure7a.png'), res=600, height=4, width=5.2, units='in')
boxplot(t(beta[,1,mc_index]), use.cols = TRUE, names=terrain_variables, horizontal=TRUE, las=1, xlab="Coefficient value", main=bquote(paste("Inflection slope, ", theta[1], sep="")))
dev.off()
png(paste0(RESULTS_DIR, '/Figure7b.png'), res=600, height=4, width=5.2, units='in')
boxplot(t(beta[,2,mc_index]), use.cols = TRUE, names=terrain_variables, horizontal=TRUE, las=1, xlab="Coefficient value", main=bquote(paste("Inflection location, ", theta[2], sep="")))
dev.off()