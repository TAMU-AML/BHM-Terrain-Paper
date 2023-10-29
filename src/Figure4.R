SOURCE_DIR = Sys.getenv("SOURCE_DIR")
LOG_DIR = Sys.getenv("LOG_DIR")
DATA_DIR = Sys.getenv("DATA_DIR")
ARTIFACTS_DIR = Sys.getenv("ARTIFACTS_DIR")
RESULTS_DIR = Sys.getenv("RESULTS_DIR")

turbine = 1
data = read.csv(paste0(DATA_DIR, '/Turbine',
                       turbine,'_2017.csv'))
thinned_ws = read.csv(paste0(DATA_DIR, '/thinned_wind_speed.csv'))[,turbine]
pdf(paste0(RESULTS_DIR, '/Figure4.pdf'), height = 4, width = 8)
par(mfrow = c(1,2))
acf(data$wind_speed, main = paste0("Turbine ", turbine," wind speed ACF"), 
    xlab = "Lag (in 10 minutes)",
    cex.lab = 1.33, cex.axis = 1.33, lag.max = 40)
acf(thinned_ws, main = paste0("Turbine ", turbine," thinned wind speed ACF"),
    xlab = "Lag (in 150 minutes)",
    cex.lab = 1.33, cex.axis = 1.33, lag.max = 40)
dev.off()
