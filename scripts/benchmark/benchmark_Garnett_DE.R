


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


# load needed libraries
library(garnett)
library(org.Mm.eg.db)


# step 1 :
# obtain DE genes from scClassify


getDE_scClassify <- function(train, trainLabels){
  
  
  doLimma <- function(exprsMat, cellTypes, exprs_pct = 0.05){
    
    cellTypes <- droplevels(as.factor(cellTypes))
    tt <- list()
    for (i in 1:nlevels(cellTypes)) {
      tmp_celltype <- (ifelse(cellTypes == levels(cellTypes)[i], 1, 0))
      design <- stats::model.matrix(~tmp_celltype)
      
      
      meanExprs <- do.call(cbind, lapply(c(0,1), function(i){
        Matrix::rowMeans(exprsMat[, tmp_celltype == i, drop = FALSE])
      }))
      
      meanPct <- do.call(cbind, lapply(c(0,1), function(i){
        Matrix::rowSums(exprsMat[, tmp_celltype == i, drop = FALSE] > 0)/sum(tmp_celltype == i)
      }))
      
      keep <- meanPct[,2] > exprs_pct
      
      y <- methods::new("EList")
      y$E <- exprsMat[keep, ]
      fit <- limma::lmFit(y, design = design)
      fit <- limma::eBayes(fit, trend = TRUE, robust = TRUE)
      
      i <- levels(cellTypes)[i]
      tt[[i]] <- limma::topTable(fit, n = Inf, adjust.method = "BH", coef = 2)
      
      
      
      if (!is.null(tt[[i]]$ID)) {
        tt[[i]] <- tt[[i]][!duplicated(tt[[i]]$ID),]
        rownames(tt[[i]]) <- tt[[i]]$ID
      }
      
      tt[[i]]$meanExprs.1 <- meanExprs[rownames(tt[[i]]), 1]
      tt[[i]]$meanExprs.2 <- meanExprs[rownames(tt[[i]]), 2]
      tt[[i]]$meanPct.1 <- meanPct[rownames(tt[[i]]), 1]
      tt[[i]]$meanPct.2 <- meanPct[rownames(tt[[i]]), 2]
    }
    
    return (tt)
    
  }
  
  tt <- doLimma(train, trainLabels)
  
  de <-  lapply(tt, function(t) rownames(t)[1:max(min(50, sum(t$adj.P.Val < 0.001)), 30)])
  
  setwd("/dora/nobackup/yuec/scclassify/benchmark")
  sink("outfile.txt")
  for (i in 1: length(de) ) {
    thislist <- de[i]
    cat(paste0(">", names(thislist)))
    cat("\n")
    cat("expressed: ")
    g = 0
    for (g in (1: (length(thislist[[1]])-1))){
      
      #print(g)
      cat(thislist[[1]][g])
      cat(", ")
    }
    cat(thislist[[1]][g+1])
    cat("\n\n")
  }
  sink()
  
}

source("/dora/nobackup/yuec/scclassify/benchmark/Garnett_DE/helper.R")

run_Garnett_DE <- function(train_data, trainLabels, test_data, testLabels, n){
  
  
  getDE_scClassify(train_data, trainLabels)
  
  
  # read the data
  cells_train <- as.matrix(colnames(train_data))
  lab_train = trainLabels
  train = as.matrix(train_data)

  
  
  # read the genefile 
  fdata <- data.frame(rownames(train_data))
  names(fdata) <- 'gene_short_name'
  row.names(fdata) <- fdata$gene_short_name
  fd_train <- new("AnnotatedDataFrame", data = fdata)
  
  
  
  true_labels <- list()
  pred_labels <- list()
  
  
  pdata_train = data.frame(cells_train)
  rownames(pdata_train) <- cells_train
  
  
  row.names(train) <- row.names(fd_train)
  colnames(train) <- row.names(pdata_train)
  
  pd_train <- new("AnnotatedDataFrame", data = pdata_train)
  pbmc_cds_train <- newCellDataSet(as(train, "dgCMatrix"), phenoData = pd_train, featureData = fd_train)
  pbmc_cds_train <- estimateSizeFactors(pbmc_cds_train)
  
  
  # training

  
  MarkerPath = "/dora/nobackup/yuec/scclassify/benchmark/outfile.txt"
  
  
  
  
  
  #use ensembl if "ENSG00000063660"
  #use SYMBOL if already converted to gene 
  marker_check <- check_markers(pbmc_cds_train, MarkerPath,
                                db=org.Mm.eg.db,
                                cds_gene_id_type = "SYMBOL",
                                marker_file_gene_id_type = "SYMBOL") #"SYMBOL"
  plot_markers(marker_check)
  
  
  benchmark <- mark(
                      pbmc_classifier <- train_cell_classifier(cds = pbmc_cds_train, 
                                           marker_file = MarkerPath,
                                           db=org.Mm.eg.db,
                                           cds_gene_id_type = "SYMBOL", #"ENSEMBL"
                                           marker_file_gene_id_type = "SYMBOL",
                                           min_observations=1,
                                           max_training_samples=50000,
                                           propogate_markers = TRUE,
                                           cores=1,
                                           lambdas = NULL,
                                           classifier_gene_id_type = "SYMBOL",
                                           return_initial_assign = FALSE),
                      time_unit = "s"
                       )
  
  
    returnlist <- list()
    returnlist$mem_train <- benchmark[, "mem_alloc"] 
    returnlist$totaltime_train <- benchmark[, "total_time"]

      
     
      lab_test = testLabels
      cells_test = as.matrix(colnames(test_data))
      test = as.matrix(test_data)
      
      # read the genefile 
      fdata <- data.frame(rownames(test_data))
      names(fdata) <- 'gene_short_name'
      row.names(fdata) <- fdata$gene_short_name
      fd_test <- new("AnnotatedDataFrame", data = fdata)
      
      
      
      true_labels <- list()
      pred_labels <- list()
      
      
      pdata_test = data.frame(cells_test)
      rownames(pdata_test) <- cells_test
      
      
      row.names(test) <- row.names(fd_test)
      colnames(test) <- row.names(pdata_test)
      
      pd_test <- new("AnnotatedDataFrame", data = pdata_test)
      pbmc_cds_test <- newCellDataSet(as(test, "dgCMatrix"), phenoData = pd_test, featureData = fd_test)
      pbmc_cds_test <- estimateSizeFactors(pbmc_cds_test)
      
      
      
      
      
      ## testing 
      
      
      benchmark <- mark (
                        pbmc_cds_test <- classify_cells(pbmc_cds_test, 
                                      pbmc_classifier, 
                                      db = org.Mm.eg.db, 
                                      cluster_extend = TRUE,
                                      cds_gene_id_type = "SYMBOL"),
                        time_unit = "s"
                          ) 
      
      returnlist$mem_test <- benchmark[, "mem_alloc"] 
      returnlist$totaltime_test <- benchmark[, "total_time"]
      
      
      true_labels <- list(lab_test)
      pred_labels <- list(pData(pbmc_cds_test)$cluster_ext_type)
      
      
      true_labels <- as.vector(unlist(true_labels))
      pred_labels <- as.vector(unlist(pred_labels))
      
      
      write.csv(true_labels, paste("train_", n, '_Garnett_DE_True.csv' , sep=""), row.names = FALSE)
      write.csv(pred_labels, paste("train_", n, '_Garnett_DE_Pred.csv', sep=""), row.names = FALSE)
      
      return (returnlist)
  }
  
  


######  Garnett_DE ########
setwd("/dora/nobackup/yuec/scclassify/benchmark/Garnett_DE")
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
  
  setwd("/dora/nobackup/yuec/scclassify/benchmark/Garnett_DE")
  returnlist <- run_Garnett_DE(train, train_label, test, test_label, thisdata )
  df[nrow(df)+1,  ] <- c(thisdata, returnlist$mem_train, returnlist$totaltime_train,
                         returnlist$mem_test, returnlist$totaltime_test)
  
}


setwd("/dora/nobackup/yuec/scclassify/benchmark/Garnett_DE")
write.csv(df, "cpu_mem.csv")


