
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



library(scmap)
library(SingleCellExperiment)


run_scmapcell <- function(train , test , train_label , test_label , n   ){
  
  True_Labels_scmapcluster <- list()
  Pred_Labels_scmapcluster <- list()
  
  
  
  trainData =  as.matrix(train)
  testData =  as.matrix(test)
  
  
  
  #prepare the training dataset
  train_label <- as.data.frame(train_label)
  colnames( train_label) <- "cell_type1"
  sce <- SingleCellExperiment(list(logcounts = trainData), 
                              colData = DataFrame(train_label))
  sce$cellTypes <- sce$cell_type1
  
  
  # use gene names as feature symbols
  rowData(sce)$feature_symbol <- rownames(sce)
  sce <- selectFeatures(sce, suppress_plot = TRUE)
  
  
  #prepare the testing dataset
  test_label <- as.data.frame(test_label)
  colnames( test_label) <- "cell_type1"
  sce_test <- SingleCellExperiment(list(logcounts = testData), 
                                   colData = DataFrame(test_label))
  sce_test$cellTypes <- sce_test$cell_type1
  
  rowData(sce_test)$feature_symbol <- rownames(sce_test)
  sce_test <- selectFeatures(sce_test, suppress_plot = TRUE)
  
  #sce_test@rowRanges@elementMetadata@listData = sce@rowRanges@elementMetadata@listData
  
  
  
  
  #scmap cell 
  set.seed(1)
  
  benchmark <- mark( sce <- indexCell(sce)  ,
                     time_unit = "s")
  
  returnlist <- list()
  returnlist$mem_train <- benchmark[, "mem_alloc"] 
  returnlist$totaltime_train <- benchmark[, "total_time"]
  
  
  
  testfunc <- function(sce_test , sce) {
    
    scmapCell_results <- scmapCell(sce_test, list(sce@metadata$scmap_cell_index))
    scmapCell_clusters <- scmapCell2Cluster(scmapCell_results,list(as.character(colData(sce)$cell_type1)))
    
    return (scmapCell_clusters)
  } 
  
  
  benchmark <- mark (scmapCell_clusters <- testfunc(sce_test, sce), 
                     time_unit = "s")
  returnlist$mem_test <- benchmark[, "mem_alloc"] 
  returnlist$totaltime_test <- benchmark[, "total_time"]
  
  
  
  
  True_Labels_scmapcell <- colData(sce_test)$cellTypes
  Pred_Labels_scmapcell <- list(scmapCell_clusters$combined_labs)
  
  True_Labels_scmapcell <- as.vector(unlist(True_Labels_scmapcell))
  Pred_Labels_scmapcell <- as.vector(unlist(Pred_Labels_scmapcell))
  
  write.csv(True_Labels_scmapcell,paste('test_' , n , '_scmapcell_True.csv',sep=""),row.names = FALSE)
  write.csv(Pred_Labels_scmapcell,paste('test_' , n , '_scmapcell_Pred.csv',sep=""),row.names = FALSE)
  
  return (returnlist)
}






######  scmapcell    ########
setwd("/dora/nobackup/yuec/scclassify/benchmark/scmapcell")
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
  
  setwd("/dora/nobackup/yuec/scclassify/benchmark/scmapcell")
  returnlist <- run_scmapcell(test, train,  test_label, train_label, thisdata )
  
  df[nrow(df)+1,  ] <- c(thisdata, returnlist$mem_train, returnlist$totaltime_train,
                         returnlist$mem_test, returnlist$totaltime_test)
  
}

setwd("/dora/nobackup/yuec/scclassify/benchmark/scmapcell")
write.csv(df, "cpu_mem.csv")




