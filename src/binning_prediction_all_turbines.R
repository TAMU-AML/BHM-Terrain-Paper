SOURCE_DIR = Sys.getenv("SOURCE_DIR")
LOG_DIR = Sys.getenv("LOG_DIR")
DATA_DIR = Sys.getenv("DATA_DIR")
ARTIFACTS_DIR = Sys.getenv("ARTIFACTS_DIR")
RESULTS_DIR = Sys.getenv("RESULTS_DIR")

source(paste0(SOURCE_DIR, '/f_binning.R'))

for (test_turb in 1:66){
  test_data = read.csv(paste0(DATA_DIR, '/Turbine',test_turb,'_2017.csv'))
  test_x = test_data$wind_speed
  test_pred = matrix(nrow=nrow(test_data), ncol=66)
  for (train_turb in 1:66){
    train_data = read.csv(paste0(DATA_DIR, '/Turbine',train_turb,'_2017.csv'))
    train_x = train_data$wind_speed
    train_y = train_data$power
    test_pred[,train_turb] = binning(train.y=train_y, train.x=train_x, test.x=test_x)
  }
  write.table(test_pred, file=paste0(ARTIFACTS_DIR, "/Turbine",test_turb,"_pred_binning.csv"), 
              sep=',', row.names=FALSE, col.names=paste0("T",c(1:66)))
}
