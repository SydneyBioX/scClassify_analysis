

.libPaths("/dora/nobackup/yuec/R")
library(bench)




####### Does scClassify’s testing time scale with training size?  #####
dataset <- c("2class", "4class", "6class", "8class", "10class", "12class")


test <- read.csv("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/test_5000.csv")
rownames(test) <- test[ , 1]
test <- test[ , -1]

test_label <- read.csv("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/test_5000_label.csv")
test_label <- as.character(test_label$x)





library(CHETAH)
library(SingleCellExperiment)
setwd("/dora/nobackup/yuec/scclassify/benchmark/CHETAH")




run_CHETAH <-function(traindata, train_label , testdata,  test_label, n ){
  
  #############################################################################
  #                                CHETAH                                     #
  #############################################################################
  
  True_Labels_CHETAH <- list()
  Pred_Labels_CHETAH <- list()
  
  
  
  sce <- SingleCellExperiment(assays = list(counts = as.matrix(traindata)), 
                              colData = DataFrame(train_label))
  
  colData(sce)$celltypes <- as.character(colData(sce)$train_label)
  #colData(sce)$ref_types <- as.character(colData(sce)$train_label)
  
  sce$celltypes
  sce_test <- SingleCellExperiment(assays = list(counts = as.matrix(testdata)),
                                   colData = DataFrame(test_label))
  
  
  colData(sce_test)$celltypes <- as.character(colData(sce_test)$test_label)
  # colData(sce_test)$ref_types <- as.character(colData(sce_test)$test_label)
  
  
  benchmark <- mark(sce_test <- CHETAHclassifier(input = sce_test, ref_cells = sce),
                    time_unit = "s")
  
  returnlist <- list()
  returnlist$mem <- benchmark[, "mem_alloc"] 
  returnlist$totaltime <- benchmark[, "total_time"]
  #https://stackoverflow.com/questions/51295402/r-on-macos-error-vector-memory-exhausted-limit-reached
  
  True_Labels_CHETAH <- colData(sce_test)$celltypes
  Pred_Labels_CHETAH <- list(sce_test$celltype_CHETAH)
  
  
  True_Labels_CHETAH <- as.vector(unlist(True_Labels_CHETAH))
  Pred_Labels_CHETAH <- as.vector(unlist(Pred_Labels_CHETAH))
  
  
  write.csv(True_Labels_CHETAH,paste("train_", n, '_CHETAH_True.csv',sep=""),row.names = FALSE)
  write.csv(Pred_Labels_CHETAH,paste("train_",  n ,  '_CHETAH_Pred.csv',sep=""),row.names = FALSE)
  
  return (returnlist)
}



df <- data.frame(n = integer(), mem = double(), totaltime = double())

for (i in (1:length(dataset))){
  thisdata <- dataset[i]
  print(thisdata)
  train <- read.csv(paste0("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/train_", thisdata, ".csv"))
  rownames(train) <- train[ , 1]
  train <- train[ , -1]
  train_label <- read.csv(paste0("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris_benchmark/train_", thisdata, "_label.csv"))
  train_label <- as.character(train_label$x)
  
  setwd("/dora/nobackup/yuec/scclassify/benchmark/CHETAH")
  tryCatch( {
  returnlist <- run_CHETAH(train, train_label, test, test_label, thisdata )
  df[nrow(df)+1,  ] <- c(thisdata, returnlist$mem, returnlist$totaltime)
  },
  error = function(e){
    
  }
  
  ) 
}

setwd("/dora/nobackup/yuec/scclassify/benchmark/CHETAH/vary_celltype")
write.csv(df, "cpu_mem.csv")