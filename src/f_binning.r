## Function to predict response from IEC binning method 
binning = function(train.x, train.y, test.x, bin_width = 0.5){
  library(stats)
  train.y = as.numeric(train.y)
  train.x = as.numeric(train.x)
  test.x = as.numeric(test.x)
  start = 0
  end = round(max(train.x))
  n_bins = round((end - start)/bin_width,0)+1
  x_bin = 0
  y_bin = 0
  for (n in 2:n_bins){
    bin_element = which(train.x>(start+(n-1)*bin_width) & train.x<(start+n*bin_width))
    x_bin[n] = mean(train.x[bin_element])
    y_bin[n] = mean(train.y[bin_element])
  }
  binned_data = data.frame(x_bin, y_bin)
  binned_data = binned_data[which(is.na(binned_data$y_bin)==F),]
  splinefit = smooth.spline(x = binned_data$x_bin, y= binned_data$y_bin , all.knots = T)
  y_pred = predict(splinefit,test.x)$y
  y_pred[which(y_pred<0)] = 0
  return(y_pred) 
}