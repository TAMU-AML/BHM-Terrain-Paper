SOURCE_DIR = Sys.getenv("SOURCE_DIR")
LOG_DIR = Sys.getenv("LOG_DIR")
DATA_DIR = Sys.getenv("DATA_DIR")
ARTIFACTS_DIR = Sys.getenv("ARTIFACTS_DIR")
RESULTS_DIR = Sys.getenv("RESULTS_DIR")

source(paste0(SOURCE_DIR, '/BHM-new-subroutines.R'))
test_turbine = 1
test_data = read.csv(paste0(DATA_DIR, '/Turbine',test_turbine,'_2017.csv'))
terrain_all = as.matrix(read.csv(paste0(DATA_DIR, "/weightedTerrainData.csv"))[c(2:4)])  
terrain_level = read.csv(paste0(DATA_DIR, '/terrain_level_agg.csv'))[,5]
turbine_locations = as.matrix(read.csv(paste0(DATA_DIR, '/location.csv'))[,-1])
terrain_cols = c(1:3)
test_wind_speed = as.matrix(test_data$wind_speed)
test_temperature = as.matrix(test_data$temperature)
test_terrain = terrain_all[test_turbine, terrain_cols, drop=FALSE]
test_power = as.matrix(test_data$power)
bhm_output = readRDS(paste0(ARTIFACTS_DIR, "/BHM_test_turbine_",test_turbine,"_all_terrain_vars.rds"))
test_pred_with_terrain = predict(bhm_output, test_wind_speed, test_temperature, test_terrain)
file_prefix = paste0(ARTIFACTS_DIR, '/Turbine')
file_suffix = '_pred_binning.csv'
file_header = TRUE
stored_pred = read.csv(paste0(file_prefix,test_turbine,file_suffix), header=file_header)
curve_types = c("Avg", "Terrain", "NN", "FN")
for (curve_type in curve_types){
  opt = list()
  opt$method = curve_type
  opt$num_neighbors = 10
  if (opt$method == "Avg"){
    train_turbines = c(1:66)[-test_turbine]
    avg_pred = rowMeans(stored_pred[ ,train_turbines, drop=FALSE])
  } else if (opt$method == "Terrain"){
    train_turbines = which(terrain_level == terrain_level[test_turbine])
    train_turbines = setdiff(train_turbines, test_turbine)
    terrain_pred = rowMeans(stored_pred[ ,train_turbines, drop=FALSE])
  } else if (opt$method == "NN"){
    neighbors = FNN::get.knnx(turbine_locations, query = turbine_locations[test_turbine,,drop=FALSE], k=nrow(turbine_locations))$nn.index[1,-1]
    train_turbines  = neighbors[1:opt$num_neighbors]
    nn10_pred = rowMeans(stored_pred[ ,train_turbines, drop=FALSE])
  } else if (opt$method == "FN"){
    neighbors = FNN::get.knnx(turbine_locations, query = turbine_locations[test_turbine,,drop=FALSE], k=nrow(turbine_locations))$nn.index[1,-1]
    train_turbines  = neighbors[(length(neighbors)-opt$num_neighbors+1):length(neighbors)]
    fn10_pred = rowMeans(stored_pred[ ,train_turbines, drop=FALSE])
  }
  
}
train_pred = stored_pred[ , test_turbine]

min_pw = min(test_power) + 20
max_pw = max(test_power) - 20
set.seed(1)
random_points = runif(100, min_pw, max_pw)
nn = FNN::get.knnx(test_power, query = random_points, k=1)
sample_idx = nn$nn.index
pdf("results/Figure9b.pdf", height = 6, width = 6)
plot(test_pred_with_terrain[sample_idx], test_power[sample_idx], cex=0.5, col="blue", xlab="Predicted power", ylab="Actual power")
points(seq(0,100,1), seq(0,100,1), type='l', lty=2, lwd=2)
points(avg_pred[sample_idx], test_power[sample_idx], pch=2, cex=0.5, col="red")
points(terrain_pred[sample_idx], test_power[sample_idx], pch=3, cex=0.5, col="green")
points(nn10_pred[sample_idx], test_power[sample_idx], pch=4, cex=0.5, col="brown")
points(fn10_pred[sample_idx], test_power[sample_idx], pch=5, cex=0.5, col="magenta")
legend("topleft", legend = c("BHM","Avg-binning","Terrain-binning","10NN-binning","10FN-binning"), pch = c(1,2,3,4,5), col =  c("blue","red","green", "brown", "magenta"), cex = 0.75)
dev.off()


sample_idx = which(test_temperature > -15 & test_temperature < -14)
pdf("results/Figure9a.pdf", height = 6, width = 6)
plot(test_wind_speed[sample_idx], test_power[sample_idx], pch=8, cex=0.5, col="black", , xlab="Wind speed (m/s)", ylab="Normalized power")
points(test_wind_speed[sample_idx], test_pred_with_terrain[sample_idx], pch=1, cex=0.5, col="blue")
points(test_wind_speed[sample_idx], avg_pred[sample_idx], pch=2, cex=0.5, col="red")
points(test_wind_speed[sample_idx], terrain_pred[sample_idx], pch=3, cex=0.5, col="green")
points(test_wind_speed[sample_idx], nn10_pred[sample_idx], pch=4, cex=0.5, col="brown")
points(test_wind_speed[sample_idx], fn10_pred[sample_idx],  pch=5, cex=0.5, col="magenta")
legend("topleft", legend = c("Data", "BHM","Avg-binning","Terrain-binning","10NN-binning","10FN-binning"), pch = c(8,1,2,3,4,5), col =  c("black","blue","red","green", "brown", "magenta"), cex = 0.75)
dev.off()


