SOURCE_DIR = Sys.getenv("SOURCE_DIR")
LOG_DIR = Sys.getenv("LOG_DIR")
DATA_DIR = Sys.getenv("DATA_DIR")
ARTIFACTS_DIR = Sys.getenv("ARTIFACTS_DIR")
RESULTS_DIR = Sys.getenv("RESULTS_DIR")

option_list = list(
  optparse::make_option("--test_turb", type="integer", default=NULL, metavar="character"),
  optparse::make_option("--terrain_vars", type="character", default="all", metavar="character")
)

opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)

test_turb = opt$test_turb
terrain_vars = opt$terrain_vars

if (terrain_vars == "all"){
  terrain_cols = c(1:3)
} else {
  terrain_cols = as.numeric(strsplit(terrain_vars, split=",")[[1]])
}

train_turbines = setdiff(c(1:66), test_turb)

terrain_all = as.matrix(read.csv(paste0(DATA_DIR,"/weightedTerrainData.csv"))[c(2:4)])  


terrain = terrain_all[train_turbines, terrain_cols, drop=FALSE] #Remove test_turbines and extra covariates

wind_speed = as.matrix(read.csv(paste0(DATA_DIR, '/thinned_wind_speed.csv')))
wind_speed = wind_speed[,train_turbines, drop=FALSE]

power = as.matrix(read.csv(paste0(DATA_DIR, '/thinned_power.csv')))
power = power[,train_turbines, drop=FALSE]

temperature = as.matrix(read.csv(paste0(DATA_DIR, '/thinned_temperature.csv')))
temperature = temperature[,train_turbines, drop=FALSE]

source(paste0(SOURCE_DIR, '/BHM-new-subroutines.R'))

start_time = Sys.time()
bhm_output = bhm_mcmc_computation(power, wind_speed, temperature, terrain, rated_power=100, burn=5000, nmc=10000, thin=30)
end_time = Sys.time()
compute_time = as.numeric(difftime(end_time, start_time, units='secs'))
cat('Time taken to compute:', compute_time, 'seconds.\n')

bhm_output$config$train_turbines = train_turbines
bhm_output$config$test_turbine = test_turb
bhm_output$config$terrain_cols = terrain_cols

saveRDS(bhm_output, file=paste0(ARTIFACTS_DIR, "/BHM_test_turbine_",test_turb,"_",terrain_vars,"_terrain_vars.rds"))  

test_data = read.csv(paste0(DATA_DIR, '/Turbine',test_turb,'_2017.csv'))
test_wind_speed = as.matrix(test_data$wind_speed)
test_temperature = as.matrix(test_data$temperature)
test_terrain = terrain_all[test_turb, terrain_cols, drop=FALSE]
test_power = as.matrix(test_data$power)
test_pred = predict(bhm_output, test_wind_speed, test_temperature, test_terrain)

rmse = sqrt(mean((test_power - test_pred)**2))
mae = mean(abs(test_power - test_pred))

bhm_output$config$test_rmse = rmse
bhm_output$config$test_mae = mae
bhm_output$config$compute_time = compute_time


saveRDS(bhm_output, file=paste0(ARTIFACTS_DIR, "/BHM_test_turbine_",test_turb, "_", terrain_vars,"_terrain_vars.rds"))  
