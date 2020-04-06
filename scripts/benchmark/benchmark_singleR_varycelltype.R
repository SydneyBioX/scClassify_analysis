

.libPaths("/dora/nobackup/yuec/R")
library(bench)




####### Does scClassify’s testing time scale with training size?  #####
dataset <- c("2class", "4class", "6class", "8class", "10class", "12class")


test <- read.csv("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/test_5000.csv")
rownames(test) <- test[ , 1]
test <- test[ , -1]

test_label <- read.csv("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/test_5000_label.csv")
test_label <- as.character(test_label$x)




library(SingleR)
library(Seurat)


run_SingleR<-function(traindata, testdata,  trainLabels, testLabels, n){
  
  
  trainLabels <- as.vector(trainLabels)
  testLabels <- as.vector(testLabels)
  
  
  
  True_Labels_SingleR <- list()
  Pred_Labels_SingleR <- list()
  
  
  
  
  traindata = as.matrix(traindata)
  testdata = as.matrix(testdata)
  
  
  
  
  Rprofmem("/dora/nobackup/yuec/scclassify/benchmark/singleR/Rprofmem.out")
  
  start.time <- Sys.time()
  singler <- SingleR(method = "single", testdata, traindata, trainLabels, numCores = 1)
  end.time <- Sys.time() 
  
  time.taken <- difftime(end.time , start.time , units = "secs")
  time.taken
  
  Rprofmem(NULL)
  mem_alloc <- sum(as.numeric(read.table("/dora/nobackup/yuec/scclassify/benchmark/singleR/Rprofmem.out", comment.char = ":", sep= ",")[,1]), na.rm=TRUE)
  
  
  
  returnlist <- list()
  returnlist$mem  <-  mem_alloc
  returnlist$totaltime  <- as.numeric(time.taken)
  
  
  
  True_Labels_SingleR <- list(testLabels)
  Pred_Labels_SingleR <- list(as.vector(singler$labels))
  
  
  True_Labels_SingleR <- as.vector(unlist(True_Labels_SingleR))
  Pred_Labels_SingleR <- as.vector(unlist(Pred_Labels_SingleR))
  
  
  
  write.csv(True_Labels_SingleR,paste("train_", n, '_SingleR_True.csv',sep=""), row.names = FALSE)
  write.csv(Pred_Labels_SingleR,paste("train_" , n ,'_SingleR_Pred.csv',sep=""), row.names = FALSE)
  
  return (returnlist)
  
}




##### SingleR  ############

setwd("/dora/nobackup/yuec/scclassify/benchmark/singleR")
df <- data.frame(n = integer(), mem = double(), totaltime = double())

for (i in (1:length(dataset))){
  thisdata <- dataset[i]
  print(thisdata)
  train <- read.csv(paste0("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/train_", thisdata, ".csv"))
  rownames(train) <- train[ , 1]
  train <- train[ , -1]
  train_label <- read.csv(paste0("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/train_", thisdata, "_label.csv"))
  train_label <- as.character(train_label$x)
  
  setwd("/dora/nobackup/yuec/scclassify/benchmark/singleR")
  tryCatch( {
    returnlist <- run_SingleR(train,  test, train_label,  test_label, thisdata )
    df[nrow(df)+1,  ] <- c(thisdata, returnlist$mem, returnlist$totaltime)
  },
  error = function(e){
    
  }
  
  ) 
  
}
setwd("/dora/nobackup/yuec/scclassify/benchmark/singleR/vary_celltype")
write.csv(df, "cpu_mem.csv")