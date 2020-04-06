


.libPaths("/dora/nobackup/yuec/R")
library(bench)
library(devtools)
install_github("satijalab/seurat", ref = "65b77a9")



####### Does scClassify’s testing time scale with training size?  #####
dataset <- c(100, 200, 500, 1000, 2000, 5000,
             10000, 20000, 30000)


test <- read.csv("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/test_2000.csv")
rownames(test) <- test[ , 1]
test <- test[ , -1]

test_label <- read.csv("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/test_2000_label.csv")
test_label <- as.character(test_label$x)



library(scPred)
library(Seurat)
library(tidyverse)
library(SingleCellExperiment)

run_scPred<-function(train , test , train_label, test_label , n ){
  
  
  True_Labels_scPred <- list()
  Pred_Labels_scPred <- list()
  
  trainData =  as.matrix( train )
  testData = as.matrix( test )  
  
  train_label <- as.data.frame(train_label)
  colnames(train_label ) <- "cellTypes"
  sce <- SingleCellExperiment(list(logcounts = trainData ), 
                              colData = DataFrame(cell_type1 = train_label))
  
  sce_counts <- logcounts(sce)
  sce_cpm <- apply(sce_counts, 2, function(x) (x/sum(x))*1000000)
  sce_metadata <- as.data.frame(colData(sce))
  
  
  # scPred Training    
  set.seed(1234)
  
  trainfunc <- function(sce_cpm,   sce_metadata){
    scp <- eigenDecompose(sce_cpm)
    scPred::metadata(scp) <- sce_metadata
    scp <- getFeatureSpace(scp, pVar = 'cell_type1.cellTypes')
    scp <- trainModel(scp)
    
    return (scp)
  }
  
  benchmark  <- mark( scp <- trainfunc (sce_cpm , sce_metadata),
                      time_unit = "s")
  
  
  returnlist <- list()
  returnlist$mem_train <- benchmark[, "mem_alloc"] 
  returnlist$totaltime_train <- benchmark[, "total_time"]
  
  
  
  
  
  
  test_label <- as.data.frame(test_label)  
  colnames(test_label) <- "cellTypes"    
  sce_test <- SingleCellExperiment(list(logcounts = testData ), 
                                   colData = DataFrame(cell_type1 = test_label))
  
  
  sce_counts_test <- logcounts(sce_test)
  sce_cpm_test <- apply(sce_counts_test, 2, function(x) (x/sum(x))*1000000)
  sce_metadata_test <- as.data.frame(colData(sce_test))
  
  
  # scPred Prediction
  
  
  benchmark <- mark (  scp <- scPredict(scp,newData = sce_cpm_test) ,
                       time_unit = "s")
  
  returnlist$mem_test <- benchmark[, "mem_alloc"] 
  returnlist$totaltime_test <- benchmark[, "total_time"]
  
  
  
  
  True_Labels_scPred <- colData(sce_test)$cell_type1.cellTypes
  Pred_Labels_scPred <- list(getPredictions(scp)$predClass)
  
  True_Labels_scPred <- as.vector(unlist(True_Labels_scPred))
  Pred_Labels_scPred <- as.vector(unlist(Pred_Labels_scPred))
  
  write.csv(True_Labels_scPred,paste("test", n, '_scPred_True.csv',sep=""), row.names = FALSE)
  write.csv(Pred_Labels_scPred,paste("test", n, '_scPred_Pred.csv',sep=""), row.names = FALSE)
  
  return (returnlist)
  
}


setwd("/dora/nobackup/yuec/scclassify/benchmark/scPred")
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
  
  setwd("/dora/nobackup/yuec/scclassify/benchmark/scPred")
  
  tryCatch( {
    
    returnlist <- run_scPred(test,  train, test_label,  train_label, thisdata )
    df[nrow(df)+1,  ] <- c(thisdata, returnlist$mem_train, returnlist$totaltime_train,
                           returnlist$mem_test, returnlist$totaltime_test)
    
  },
  error = function(e){
  }
  
  )
  
}

setwd("/dora/nobackup/yuec/scclassify/benchmark/scPred")
write.csv(df, "cpu_mem.csv")


