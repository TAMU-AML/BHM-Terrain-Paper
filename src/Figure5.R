SOURCE_DIR = Sys.getenv("SOURCE_DIR")
LOG_DIR = Sys.getenv("LOG_DIR")
DATA_DIR = Sys.getenv("DATA_DIR")
ARTIFACTS_DIR = Sys.getenv("ARTIFACTS_DIR")
RESULTS_DIR = Sys.getenv("RESULTS_DIR")

powercurve = function(x, theta){
  return(1/(1+exp(-theta[1]*(x-theta[2]))))
}

x = seq(0, 20, length.out = 100)

theta1 = c(.7, 8)
f1 = powercurve(x, theta1)

theta2 = c(1.2, 8)
f2 = powercurve(x, theta2)

theta3 = c(.7, 8)
f3 = powercurve(x, theta3)

theta4 = c(.7, 9)
f4 = powercurve(x, theta4)

pdf(paste0(RESULTS_DIR, "/Figure5.pdf"), height = 4, width = 8)
par(mfrow=c(1,2))
plot(x, f1, type = 'l', lty = 1, lwd = 2, col = "blue",
     main = expression("a) Increasing"~theta[1]), 
     xlab = "Wind speed (m/s)",
     ylab = "Normalized power")
points(x, f2, type = 'l', lty = 2, lwd = 2, col = "red")

plot(x, f3, type = 'l', lty = 1, lwd = 2, col = "blue",
     main = expression("b) Increasing"~theta[2]), 
     xlab = "Wind speed (m/s)",
     ylab = "Normalized power")
points(x, f4, type = 'l', lty = 2, lwd = 2, col = "red")
dev.off()