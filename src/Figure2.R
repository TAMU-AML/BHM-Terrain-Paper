SOURCE_DIR = Sys.getenv("SOURCE_DIR")
LOG_DIR = Sys.getenv("LOG_DIR")
DATA_DIR = Sys.getenv("DATA_DIR")
ARTIFACTS_DIR = Sys.getenv("ARTIFACTS_DIR")
RESULTS_DIR = Sys.getenv("RESULTS_DIR")

locations = read.csv(paste0(DATA_DIR, '/location.csv'))
terrain_level = read.csv(paste0(DATA_DIR, '/terrain_level_agg.csv'))
{
  pch = rep(NA, 66)
  pch[which(terrain_level[,5] == 3)] = 19
  pch[which(terrain_level[,5] == 4)] = 17
  pch[which(terrain_level[,5] == 5)] = 15
  
  
  png(paste0(RESULTS_DIR, '/Figure2.png'), res=600, height=4, width=8, units='in')
  par(mai=c(0.25,0.25,0.5,0.25))
  plot(locations$Longitude, locations$Latitude, pch=pch,
       main="Wind Farm Layout", cex.main=1.5, xaxt='n', yaxt='n',xlab=NA, ylab=NA)
  text(locations$Longitude, locations$Latitude, labels=c(1:66), pos=3, offset=.3, cex=0.5)
  legend("bottomright", c("Level 3 (44)", "Level 4 (5)", "Level 5 (17)"), pch=c(19, 17, 15), cex=1.33)
  dev.off()
}