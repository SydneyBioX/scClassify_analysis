


.libPaths("/dora/nobackup/yuec/R")
library(bench)




####### Does scClassify’s testing time scale with training size?  #####
dataset <- c(100, 200, 500, 1000, 2000, 5000,
             10000, 20000, 30000)


test <- read.csv("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/test_2000.csv")
rownames(test) <- test[ , 1]
test <- test[ , -1]

test_label <- read.csv("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/test_2000_label.csv")
test_label <- as.character(test_label$x)



library(igraph)
library(xgboost)





run_CaSTLe<-function(train, test, trainlabels, testlabels, n){
  returnlist <- list()
  returnlist$mem_train <- list()
  returnlist$totaltime_train <- list()
  
  returnlist$mem_test <- list()
  returnlist$totaltime_test <- list()
  
  
  
  
  True_Labels_Castle <- list()
  Pred_Labels_Castle <- list()
  
  
  BREAKS=c(-1, 0, 1, 6, Inf)
  nFeatures = 100
  
  train = t(train)
  test  = t(test)
  ds1 = train
  ds2 = test
  sourceCellTypes = as.factor(trainlabels[1:length(trainlabels)])
  targetCellTypes = as.factor(testlabels[1:length(testlabels)])
  
  
  
  # 2. Unify sets, excluding low expressed genes
  source_n_cells_counts = apply(ds1, 2, function(x) { sum(x > 0) } )
  target_n_cells_counts = apply(ds2, 2, function(x) { sum(x > 0) } )
  common_genes = intersect( colnames(ds1)[source_n_cells_counts>10], 
                            colnames(ds2)[target_n_cells_counts>10])
  remove(source_n_cells_counts, target_n_cells_counts)
  ds1 = ds1[, colnames(ds1) %in% common_genes]
  ds2 = ds2[, colnames(ds2) %in% common_genes]
  ds = rbind(ds1[,common_genes], ds2[,common_genes])
  isSource = c(rep(TRUE,nrow(ds1)), rep(FALSE,nrow(ds2)))
  remove(ds1, ds2)
  
  # 3. Highest mean in both source and target
  topFeaturesAvg = colnames(ds)[order(apply(ds, 2, mean), decreasing = T)]
  
  # for each cell - what is the most probable classification?
  L = length(levels(sourceCellTypes))
  targetClassification = as.data.frame(matrix(rep(0,L*sum(!isSource)), nrow=L), row.names = levels(sourceCellTypes))
  
  for (cellType in levels(sourceCellTypes)) {
    
    inSourceCellType = as.factor(ifelse(sourceCellTypes == cellType, cellType, paste0("NOT",cellType)))
    
    # 4. Highest mutual information in source
    topFeaturesMi = names(sort(apply(ds[isSource,],2,function(x) { compare(cut(x,breaks=BREAKS),inSourceCellType,method = "nmi") }), decreasing = T))
    
    # 5. Top n genes that appear in both mi and avg
    selectedFeatures = union(head(topFeaturesAvg, nFeatures) , head(topFeaturesMi, nFeatures) )
    
    # 6. remove correlated features
    tmp = cor(ds[,selectedFeatures], method = "pearson")
    tmp[!lower.tri(tmp)] = 0
    selectedFeatures = selectedFeatures[apply(tmp,2,function(x) any(x < 0.9))]
    remove(tmp)
    
    # 7,8. Convert data from continous to binned dummy vars
    # break datasets to bins
    dsBins = apply(ds[, selectedFeatures], 2, cut, breaks= BREAKS)
    # use only bins with more than one value
    nUniq = apply(dsBins, 2, function(x) { length(unique(x)) })
    # convert to dummy vars
    ds0 = model.matrix(~ . , as.data.frame(dsBins[,nUniq>1]))
    remove(dsBins, nUniq)
    
    cat(paste0("<h2>Classifier for ",cellType,"</h2>"))
    
    inTypeSource = sourceCellTypes == cellType
    # 9. Classify
    
    
    
    
    benchmark <- mark( xg <- xgboost(data=ds0[isSource,] , 
                                     label=inTypeSource,
                                     objective="binary:logistic", 
                                     eta=0.7 , nthread=1, nround=20, verbose=0,
                                     gamma=0.001, max_depth=5, min_child_weight=10),
                       time_unit = "s"
                    )
    
    
    returnlist$mem_train <- c(returnlist$mem_train,  benchmark[, "mem_alloc"] )
    returnlist$totaltime_train <- c(returnlist$totaltime_train,  benchmark[, "total_time"])
    
    
    # 10. Predict
    benchmark <- mark ( inTypeProb <- predict(xg, ds0[!isSource, ]) ,
                        time_unit = "s"
                       )
    
    returnlist$mem_test <- c(returnlist$mem_test,  benchmark[, "mem_alloc"] )
    returnlist$totaltime_test <- c(returnlist$totaltime_test,  benchmark[, "total_time"])
    
    
    
    targetClassification[cellType,] = inTypeProb
  }
  
  
  
  True_Labels_Castle <- list(testlabels[1:length(testlabels)])
  Pred_Labels_Castle <- list(rownames(targetClassification)[apply(targetClassification,2,which.max)])
  
  True_Labels_Castle <- as.vector(unlist(True_Labels_Castle))
  Pred_Labels_Castle <- as.vector(unlist(Pred_Labels_Castle))
  
  
  write.csv(True_Labels_Castle,paste("train_", n, '_True_CaSTLe.csv', sep= ""),row.names = FALSE)
  write.csv(Pred_Labels_Castle,paste("train_", n,  '_Pred_CaSTLe.csv',sep=""), row.names = FALSE)
  
  return (returnlist)
}





######  CaSTLe ########
setwd("/dora/nobackup/yuec/scclassify/benchmark/CaSTLe")
df <- data.frame(n = integer(), mem_train = double(), totaltime_train = double(),
                 mem_test = double(), totaltime_test = double())


for (i in (1:length(dataset))){
  thisdata <- dataset[i]
  print(thisdata)
  train <- read.csv(paste0("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/train_", thisdata, ".csv"))
  rownames(train) <- train[ , 1]
  train <- train[ , -1]
  train_label <- read.csv(paste0("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/train_", thisdata, "_label.csv"))
  train_label <- as.character(train_label$x)
  
  setwd("/dora/nobackup/yuec/scclassify/benchmark/CaSTLe")
  returnlist <- run_CaSTLe(train, test,  train_label, test_label, thisdata )
  returnlist$mem_train <- sum(unlist(returnlist$mem_train))
  returnlist$totaltime_train <- sum(unlist(returnlist$totaltime_train))
  returnlist$mem_test <- sum(unlist(returnlist$mem_test))
  returnlist$totaltime_test <- sum(unlist(returnlist$totaltime_test))
  
  
  df[nrow(df)+1,  ] <- c(thisdata, returnlist$mem_train, returnlist$totaltime_train,
                         returnlist$mem_test, returnlist$totaltime_test)
  
}

setwd("/dora/nobackup/yuec/scclassify/benchmark/CaSTLe")
write.csv(df, "cpu_mem.csv")


