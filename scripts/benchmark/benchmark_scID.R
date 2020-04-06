
.libPaths("/dora/nobackup/yuec/R")
library(bench)

#install.packages("Seurat")
#install.packages("MAST")


####### Does scClassify’s testing time scale with training size?  #####
dataset <- c(100, 200, 500, 1000, 2000, 5000,
             10000, 20000, 30000)


test <- read.csv("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/test_2000.csv")
rownames(test) <- test[ , 1]
test <- test[ , -1]

test_label <- read.csv("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/test_2000_label.csv")
test_label <- as.character(test_label$x)





library(scID)
library(Seurat)


run_scID<-function(train, test, trainLabels, testLabels,  n){

  True_Labels_scID <- list()
  Pred_Labels_scID <- list()

  
  trainData = train
  testData = test
  
  Train_Labels <- list(trainLabels)
  names(Train_Labels[[1]]) <- colnames(trainData) #do some formatting, so that it can be used to run scID 
  
  
  
  
  
  benchmark <- mark( scID_output <- scid_multiclass(testData, trainData, Train_Labels[[1]]) ,
                     time_unit = "s")
  
  returnlist <- list()
  returnlist$mem <- benchmark[, "mem_alloc"] 
  returnlist$totaltime <- benchmark[, "total_time"]
  
  
  True_Labels_scID <- list(testLabels)
  Pred_Labels_scID <- list(as.vector(scID_output$labels))
  True_Labels_scID <- as.vector(unlist(True_Labels_scID))
  Pred_Labels_scID <- as.vector(unlist(Pred_Labels_scID))
  
  
  write.csv(True_Labels_scID,paste("train_", n, "_scID_True.csv",sep=""), row.names = FALSE)
  write.csv(Pred_Labels_scID,paste("train_", n, "_scID_Pred.csv",sep=""), row.names = FALSE)
  
  return (returnlist)
}




###### scID ######

df <- data.frame(n = integer(), mem = double(), totaltime = double())

for (i in (1:length(dataset))){
  thisdata <- dataset[i]
  print(thisdata)
  train <- read.csv(paste0("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/train_", thisdata, ".csv"))
  rownames(train) <- train[ , 1]
  train <- train[ , -1]
  train_label <- read.csv(paste0("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/train_", thisdata, "_label.csv"))
  train_label <- as.character(train_label$x)
  
  setwd("/dora/nobackup/yuec/scclassify/benchmark/scID")
  returnlist <- run_scID(train,  test, train_label, test_label, thisdata )
  df[nrow(df)+1,  ] <- c(thisdata, returnlist$mem, returnlist$totaltime)
}

setwd("/dora/nobackup/yuec/scclassify/benchmark/scID")
write.csv(df, "cpu_mem.csv")

