

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


library(singleCellNet)
library(dplyr)

run_singleCellNet<-function(trainData, testData, trainLabels, testLabels, n){
  
  
  trainLabels <- as.vector(trainLabels)
  testLabels <- as.vector(testLabels)
  
  
  True_Labels_singleCellNet <- list()
  Pred_Labels_singleCellNet <- list()
  
  
  trainData =  as.matrix(trainData)
  testData =  as.matrix(testData)
  
  
  
  trainfunc <- function(trainData, trainLabels) {
    cgenes2  <-findClassyGenes( trainData, data.frame(Annotation = trainLabels[1:length(trainLabels)]), "Annotation")
    cgenesA  <-cgenes2[['cgenes']]
    grps  <-cgenes2[['grps']]
    
    
    trainData <-as.matrix(trainData[cgenesA,])
    xpairs<-ptGetTop(trainData, grps, ncores = 1)
    
    pdTrain<-query_transform(trainData[cgenesA, ], xpairs)
    rf <-sc_makeClassifier(pdTrain[xpairs,], genes=xpairs, groups=grps)
    
    result  <- list()
    result$cgenesA <- cgenesA
    result$rf <- rf
    result$xpairs <- xpairs
    
    return (result)
  }
  
  benchmark <- mark ( result <- trainfunc(trainData, trainLabels) ,
                      time_unit = "s" ) 
  returnlist <- list()
  returnlist$mem_train <- benchmark[, "mem_alloc"] 
  returnlist$totaltime_train <- benchmark[, "total_time"]
  
 
  
  
  testfunc <- function(testData, rf, xpairs, cgenesA ){
    testData <- testData[rownames(testData) %in% cgenesA, ]
    testData <- query_transform(testData, xpairs)
    classRes <-rf_classPredict(rf, testData)
    
    return (classRes)
  }
 
  
  benchmark <- mark (classRes  <- testfunc( testData,  result$rf, result$xpairs , result$cgenesA  ) , 
                     time_unit = "s" )
  
  returnlist$mem_test <- benchmark[, "mem_alloc"] 
  returnlist$totaltime_test <- benchmark[, "total_time"]
  
  
  
  
  True_Labels_singleCellNet <- list(testLabels)
  Pred_Labels_singleCellNet <- list((rownames(classRes)[apply(classRes,2,which.max)]))

  
  True_Labels_singleCellNet <- as.vector(unlist(True_Labels_singleCellNet))
  Pred_Labels_singleCellNet <- as.vector(unlist(Pred_Labels_singleCellNet))
  
  
  write.csv(True_Labels_singleCellNet,paste( "train_", n , '_singleCellNet_True.csv', sep=""),row.names = FALSE)
  write.csv(Pred_Labels_singleCellNet,paste("train_", n ,  '_singleCellNet_Pred.csv',sep= ""),row.names = FALSE)
  
  return (returnlist)
}




##### SingleCellNet  ############

setwd("/dora/nobackup/yuec/scclassify/benchmark/singlecellnet")
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
  
  setwd("/dora/nobackup/yuec/scclassify/benchmark/singlecellnet")
  returnlist <- run_singleCellNet(train,  test, train_label,  test_label, thisdata )
  df[nrow(df)+1,  ] <- c(thisdata, returnlist$mem_train, returnlist$totaltime_train,
                         returnlist$mem_test, returnlist$totaltime_test)
  
  
}
setwd("/dora/nobackup/yuec/scclassify/benchmark/singlecellnet")
write.csv(df, "cpu_mem.csv")



