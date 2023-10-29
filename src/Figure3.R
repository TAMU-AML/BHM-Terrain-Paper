SOURCE_DIR = Sys.getenv("SOURCE_DIR")
LOG_DIR = Sys.getenv("LOG_DIR")
DATA_DIR = Sys.getenv("DATA_DIR")
ARTIFACTS_DIR = Sys.getenv("ARTIFACTS_DIR")
RESULTS_DIR = Sys.getenv("RESULTS_DIR")

data = read.csv(paste0(DATA_DIR, '/Turbine1_2017.csv'))
data$wind_direction = data$wind_direction-5
pdf(paste0(RESULTS_DIR, "/Figure3.pdf"), height = 4, width = 6)
hist(data$wind_direction, breaks = seq(-5,355,10), probability = TRUE,
     main = "Turbine 1", ylab = "Normalized frequency", xlab = "Wind direction")
dev.off()
