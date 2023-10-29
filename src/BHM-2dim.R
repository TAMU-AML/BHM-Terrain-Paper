SOURCE_DIR = Sys.getenv("SOURCE_DIR")
LOG_DIR = Sys.getenv("LOG_DIR")
DATA_DIR = Sys.getenv("DATA_DIR")
ARTIFACTS_DIR = Sys.getenv("ARTIFACTS_DIR")
RESULTS_DIR = Sys.getenv("RESULTS_DIR")

args = commandArgs(trailingOnly = TRUE)
test_turb = args[1]
if (args[2] == "all"){
  terrain_covs = c(1:3)
} else {
  terrain_covs = as.numeric(strsplit(args[2], split=",")[[1]])
}

print(paste0("Start time: ", Sys.time()))
start_time = Sys.time()

# source necessary functions
source(paste0(SOURCE_DIR, '/BHM-new-subroutines.R'))

# Import dataset
if (test_turb == "None"){
  test_turbine = NULL
} else {
  test_turbine = as.integer(test_turb)
}

train_turbines = setdiff(c(1:66), test_turbine)

terrain = as.matrix(read.csv(paste0(DATA_DIR, "/weightedTerrainData.csv"))[-1])
terrain = terrain[train_turbines, terrain_covs, drop=FALSE] #Remove test_turbines and extra covariates

wind_speed = as.matrix(read.csv(paste0(DATA_DIR, '/thinned_wind_speed.csv')))
wind_speed = wind_speed[,train_turbines]

power = as.matrix(read.csv(paste0(DATA_DIR, '/thinned_power.csv')))
power = power[,train_turbines]

temperature = as.matrix(read.csv(paste0(DATA_DIR, '/thinned_temperature.csv')))
temperature = temperature[,train_turbines]

bhm_output = bhm_mcmc_computation(power, wind_speed, temperature, terrain, rated_power=100, burn=5000, nmc=10000, thin=30)

end_time = Sys.time()
print(end_time - start_time)

saveRDS(bhm_output, file=paste0(ARTIFACTS_DIR, "/BHM-2dim-",args[1],"-",args[2],".rds"))



