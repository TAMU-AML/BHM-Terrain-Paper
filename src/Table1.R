SOURCE_DIR = Sys.getenv("SOURCE_DIR")
LOG_DIR = Sys.getenv("LOG_DIR")
DATA_DIR = Sys.getenv("DATA_DIR")
ARTIFACTS_DIR = Sys.getenv("ARTIFACTS_DIR")
RESULTS_DIR = Sys.getenv("RESULTS_DIR")

option_list = list(
  optparse::make_option("--method", type="character", metavar="character")
)
opt_parser = optparse::OptionParser(option_list=option_list)
opt = optparse::parse_args(opt_parser)
method = opt$method
##function to create nFold CV dataset list
createCVdataset = function(dataset, nFolds){
  set.seed(1)
  cv_sample = sample(nrow(dataset))
  folds = cut(seq(1:length(cv_sample)),breaks = nFolds, labels = F)
  CVdataset = vector("list", nFolds)
  for (f in 1:nFolds){
    testindex = cv_sample[which(folds == f,arr.ind = T)]
    CVdataset[[f]]$train = dataset[-testindex,]
    CVdataset[[f]]$test = dataset[testindex,]
  }
  return(CVdataset)
}


##Subset Selection as well as in-temporal average RMSE.
forwardSubsetSelection = function(CVdataset, ycol, covariates, cirCov, bestSubset = NULL, bestRMSE = Inf, method="amk", rmselist = list(), subsetlist = list(), rmseSdList = list(), verbose = FALSE){
  addCol = setdiff(covariates, bestSubset)
  ncov = length(addCol)
  if (ncov>0){
    bestCol = NULL
    rmseCov = rep(0, ncov)
    rmseCovSd = rep(0, ncov)
    for (i in 1:ncov){
      if (verbose){
        cat("Trying covariate:",covariate_names[addCol[i]],'\n', file = logfile, append = TRUE)
      }
      rmseVec = rep(0,length(CVdataset))
      covSubset = c(bestSubset,addCol[i])
      for (f in 1:length(CVdataset)){
        if (verbose){
          cat("Fold = ",f,'\n', file = logfile, append = TRUE)
        }
        train = CVdataset[[f]]$train
        test = CVdataset[[f]]$test
        if (method == "knn"){
          knnMdl = DSWE::KnnPCFit(data = train, xCol = covSubset, yCol = ycol)
          predY = DSWE::KnnPredict(knnMdl,test)  
        } else if (method == "amk"){
          if (addCol[i] == cirCov){
            cirIdx = which(covSubset == cirCov)
          } else {cirIdx = NA}
          if (length(covSubset) <= 3){
            nMultiCov = "all"
          } else {nMultiCov = 3}
          predY = DSWE::AMK(train[,covSubset, drop=FALSE], train[, ycol], test[,covSubset, drop=FALSE], cirCov = cirIdx, nMultiCov = nMultiCov, bw = 'dpi')
        } else if (method == "svm"){
          predY = DSWE::SvmPCFit(train[,covSubset, drop=FALSE], train[, ycol], test[,covSubset, drop=FALSE], kernel = "radial")
        } else if (method == "spline"){
          predY = DSWE::SplinePCFit(train, covSubset, yCol, test)
        } 
        rmseVec[f] = sqrt(mean((predY - test[,ycol])^2))/rated_power
        if (verbose){
          cat('RMSE for fold',f,'is',rmseVec[f],'\n', file = logfile, append = TRUE)
        }  
      }
      rmseCov[i] = mean(rmseVec)
      rmseCovSd[i] = sd(rmseVec)
      if (verbose){
        cat('Average RMSE using covariate',covariate_names[addCol[i]],'is: ',rmseCov[i],'(',rmseCovSd[i],')\n', file = logfile, append = TRUE)
      }
      if (rmseCov[i] < bestRMSE){
        bestRMSE = rmseCov[i]
        bestCol = addCol[i]
      }
    }
    if (length(bestCol)==0){
      bestCol = addCol[which.min(rmseCov)]
    } 
    bestSubset = c(bestSubset,bestCol)
    subsetlist = append(subsetlist, list(bestSubset))
    rmselist = append(rmselist, min(rmseCov))
    rmseSdList = append(rmseSdList, rmseCovSd[which.min(rmseCov)])
    if (verbose){
      cat("Current best subset:",covariate_names[bestSubset],'\n', file = logfile, append = TRUE)
      cat("Current best RMSE:",round(bestRMSE,4),'\n', file = logfile, append = TRUE)
    }
    retList = forwardSubsetSelection(CVdataset, ycol, covariates, cirCov, bestSubset, bestRMSE, method=method, rmselist = rmselist, subsetlist=subsetlist, rmseSdList=rmseSdList, verbose=verbose)
  } else {
    retList = list(bestRMSE = bestRMSE, bestSubset = bestSubset, rmselist = rmselist, subsetlist=subsetlist, rmseSdList=rmseSdList)
  }
  return(retList)
}

for (turbine in 1:1){
  data = read.csv(paste0(DATA_DIR, "/Turbine",turbine,"_2017.csv"))
  xCol = c(3:7)
  covariate_names = c("id", "ts", "V", "D", "T", "I", "sdD", "y")
  yCol = 8
  cirCov = 4
  rated_power = 100
  CVdataset = createCVdataset(data, 5)
  logfile = ""
  output = forwardSubsetSelection(CVdataset, yCol, xCol, cirCov, method=method, verbose = TRUE)
  saveRDS(output, file=paste0(ARTIFACTS_DIR,"/",method,"_input_subset_turbine_",turbine,".rds"))
}

result = matrix(nrow = length(xCol), ncol = 2)
row_name = c()
for (i in 1:nrow(result)) {
  row_name = c(row_name, paste(covariate_names[output$subsetlist[[i]]], collapse = ' '))
  result[i, 1] = output$rmselist[[i]]
  result[i, 2] = output$rmseSdList[[i]]
}

write.table(round(result, 4), file = paste0(RESULTS_DIR, "/Table1.txt"), sep='\t',
            row.names = row_name, col.names = c("RMSE", "Std. Err"))