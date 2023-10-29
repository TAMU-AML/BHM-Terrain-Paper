SOURCE_DIR = Sys.getenv("SOURCE_DIR")
LOG_DIR = Sys.getenv("LOG_DIR")
DATA_DIR = Sys.getenv("DATA_DIR")
ARTIFACTS_DIR = Sys.getenv("ARTIFACTS_DIR")
RESULTS_DIR = Sys.getenv("RESULTS_DIR")

methods = c("BHM", "Avg-binning", "Terrain-binning", "10NN-binning", "20NN-binning", "10FN-binning", "20FN-binning")

file_prefix = paste0(ARTIFACTS_DIR, '/Turbine')
file_suffix = '_pred_binning.csv'
file_header = TRUE
terrain_vars = "all"

terrain_level = read.csv(paste0(DATA_DIR, '/terrain_level_agg.csv'))[,5]
turbine_locations = as.matrix(read.csv(paste0(DATA_DIR, '/location.csv'))[,-1])

output = array(0, dim=c(66, length(methods), 2))

for (test_turb in 1:66){
    test_data = read.csv(paste0(DATA_DIR, '/Turbine',test_turb,'_2017.csv'))
    test_y = test_data$power
    stored_pred = read.csv(paste0(file_prefix,test_turb,file_suffix), header=file_header)
    for (method_counter in 1:length(methods)){
      if (methods[method_counter] == "BHM") {
        bhm_model = readRDS(paste0(ARTIFACTS_DIR, '/BHM_test_turbine_',test_turb,"_",terrain_vars,"_terrain_vars.rds"))  
        output[test_turb, method_counter, 1] = bhm_model$config$test_rmse
        output[test_turb, method_counter, 2] = bhm_model$config$test_mae
        next
      } else if (methods[method_counter] == "Avg-binning"){
        train_turbines = c(1:66)[-test_turb]
      } else if (methods[method_counter] == "Terrain-binning"){
        train_turbines = which(terrain_level == terrain_level[test_turb])
        train_turbines = setdiff(train_turbines, test_turb)
      } else if (methods[method_counter] == "10NN-binning"){
        num_neighbors = 10
        neighbors = FNN::get.knnx(turbine_locations, query = turbine_locations[test_turb,,drop=FALSE], k=nrow(turbine_locations))$nn.index[1,-1]
        train_turbines  = neighbors[1:num_neighbors]
      } else if (methods[method_counter] == "20NN-binning"){
        num_neighbors = 20
        neighbors = FNN::get.knnx(turbine_locations, query = turbine_locations[test_turb,,drop=FALSE], k=nrow(turbine_locations))$nn.index[1,-1]
        train_turbines  = neighbors[1:num_neighbors]
      } else if (methods[method_counter] == "10FN-binning"){
        num_neighbors = 10
        neighbors = FNN::get.knnx(turbine_locations, query = turbine_locations[test_turb,,drop=FALSE], k=nrow(turbine_locations))$nn.index[1,-1]
        train_turbines  = neighbors[(length(neighbors)-num_neighbors+1):length(neighbors)]
      } else if (methods[method_counter] == "20FN-binning"){
        num_neighbors = 20
        neighbors = FNN::get.knnx(turbine_locations, query = turbine_locations[test_turb,,drop=FALSE], k=nrow(turbine_locations))$nn.index[1,-1]
        train_turbines  = neighbors[(length(neighbors)-num_neighbors+1):length(neighbors)]
      }
      pred = rowMeans(stored_pred[ ,train_turbines, drop=FALSE])
      rmse = sqrt(mean((test_y - pred)**2))
      mae = mean(abs(test_y - pred))
      output[test_turb, method_counter, 1] = rmse
      output[test_turb, method_counter, 2] = mae
    }
}

output = apply(output, c(2,3), mean)
write.table(round(output, 3), file = paste0(RESULTS_DIR, "/Table4.txt"), sep = '\t', 
            row.names = methods, col.names = c("RMSE", "MAE"))




