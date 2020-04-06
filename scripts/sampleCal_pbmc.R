rm(list = ls())


# Package

library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(pheatmap)
library(ggthemes)

library(limma)
library(hopach)
library(ggraph)
library(igraph)
library(mixtools)
library(pbmcapply)
library(proxy)
library(minpack.lm)



library(scdney)
library(HDF5Array)

levinPBMCbenchmark <- readRDS("data/levinPBMCbenchmark.rds")

table(levinPBMCbenchmark$Method)

library(scMerge)
data("segList_ensemblGeneID")

saveRDSwithDate <- function(fileName, ...) {
  today <- gsub("-","", Sys.Date())
  fileName <- paste(paste(fileName, today, sep = "_"), ".rds", sep = "")
  saveRDS(..., file = fileName)
}



################################################
########      Smart-seq2         ###############
################################################

levinPBMCbenchmark_smartSeq <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "Smart-seq2"]

levinPBMCbenchmark_smartSeq$batch <- droplevels(levinPBMCbenchmark_smartSeq$Experiment)

logcounts(levinPBMCbenchmark_smartSeq) <- as.matrix(logcounts(levinPBMCbenchmark_smartSeq))
counts(levinPBMCbenchmark_smartSeq) <- as.matrix(counts(levinPBMCbenchmark_smartSeq))

levinPBMCbenchmark_smartSeq <- levinPBMCbenchmark_smartSeq[rowSums(counts(levinPBMCbenchmark_smartSeq))!=0,]


colsums_smartseq <- Matrix::colSums(counts(levinPBMCbenchmark_smartSeq))
nGene_smartseq <- Matrix::colSums(counts(levinPBMCbenchmark_smartSeq)!=0)

cutoff_low <- median(nGene_smartseq) - 3*mad(nGene_smartseq)

sum(nGene_smartseq < cutoff_low)



sum(nGene_smartseq > 6000)

keep <- nGene_smartseq < 6000 & colsums_smartseq > 100000
sum(keep)


levinPBMCbenchmark_smartSeq$keep_Cell <- keep

levinPBMCbenchmark_smartSeq <- runTSNE(levinPBMCbenchmark_smartSeq)

plotTSNE(levinPBMCbenchmark_smartSeq, colour_by = "keep_Cell")
plotTSNE(levinPBMCbenchmark_smartSeq, colour_by = "Experiment")


levinPBMCbenchmark_smartSeq <- levinPBMCbenchmark_smartSeq[,keep]

levinPBMCbenchmark_smartSeq <- levinPBMCbenchmark_smartSeq[rowSums(counts(levinPBMCbenchmark_smartSeq))!=0, ]

system.time(levinPBMCbenchmark_smartSeq <- scMerge(levinPBMCbenchmark_smartSeq,
                                                   # ctl = segList_ensemblGeneID$human$human_scSEG,
                                                   ctl = rownames(levinPBMCbenchmark_smartSeq),
                                                   cell_type = levinPBMCbenchmark_smartSeq$CellType,
                                                   # cell_type_match = T,
                                                   exprs = "logcounts",
                                                   kmeansK = rep(length(unique(levinPBMCbenchmark_smartSeq$CellType)),length(unique(levinPBMCbenchmark_smartSeq$batch))),
                                                   assay_name = "scMerge",
                                                   fast_svd = T))

levinPBMCbenchmark_smartSeq <- runTSNE(levinPBMCbenchmark_smartSeq)
plotTSNE(levinPBMCbenchmark_smartSeq, colour_by = "Experiment")

levinPBMCbenchmark_smartSeq <- runTSNE(levinPBMCbenchmark_smartSeq, exprs_values = "scMerge")
plotTSNE(levinPBMCbenchmark_smartSeq, colour_by = "Experiment")
plotTSNE(levinPBMCbenchmark_smartSeq, colour_by = "CellType")
# subsetIdx <- cellTypes_pbmc10k%in%c("CD4+ T Helper2", "CD4+/CD25 T Reg", "CD4+/CD45RA+/CD25- Naive T", "CD4+/CD45RO+ Memory", "CD8+ Cytotoxic T", "CD8+/CD45RA+ Naive Cytotoxic")




levinPBMCbenchmark_smartSeq$CellType2 <- as.character(droplevels(levinPBMCbenchmark_smartSeq$CellType))
levinPBMCbenchmark_smartSeq$CellType2[grep("monocyte", levinPBMCbenchmark_smartSeq$CellType2)] <- "monocyte"
levinPBMCbenchmark_smartSeq$CellType2[grep("T cell", levinPBMCbenchmark_smartSeq$CellType2)] <- "T cell"

#################################################################################
######################### level 2 ###############################################
#################################################################################


exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_smartSeq, "scMerge"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_smartSeq$CellType)
table(cellTypes_subset)


set.seed(2019)

system.time(res_sub <- mclapply(1:10, function(x) {
  
  l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), subset_test = F,  n = 20,  balance = F)
  table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)
}, mc.cores = 10))
gc()




library(pbmcapply)

n_list <- seq(20, 400, 20)

res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = F,  n = n_list[[i]], balance = F, num_test = 100)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 50))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_smartseq_level2", res_sub)


#################################################################################
######################### level 1 ###############################################
#################################################################################


exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_smartSeq, "scMerge"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- (levinPBMCbenchmark_smartSeq$CellType2)
table(cellTypes_subset)


set.seed(2019)

system.time(res_sub <- mclapply(1:10, function(x) {
  
  l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), subset_test = F,  n = 20,  balance = F)
  table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)
}, mc.cores = 10))
gc()




library(pbmcapply)

n_list <- seq(20, 400, 20)

res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = F,  n = n_list[[i]], balance = F, num_test = 100)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 50))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_smartseq_level1", res_sub)


rm("levinPBMCbenchmark_smartSeq")
gc(reset = T)

################################################
########      Cel-seq2           ###############
################################################

levinPBMCbenchmark_celSeq <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "CEL-Seq2"]

levinPBMCbenchmark_celSeq$batch <- droplevels(levinPBMCbenchmark_celSeq$Experiment)

logcounts(levinPBMCbenchmark_celSeq) <- as.matrix(logcounts(levinPBMCbenchmark_celSeq))
counts(levinPBMCbenchmark_celSeq) <- as.matrix(counts(levinPBMCbenchmark_celSeq))

levinPBMCbenchmark_celSeq <- levinPBMCbenchmark_celSeq[rowSums(counts(levinPBMCbenchmark_celSeq))!=0,]
levinPBMCbenchmark_celSeq 

system.time(levinPBMCbenchmark_celSeq <- scMerge(levinPBMCbenchmark_celSeq,
                                                 # ctl = segList_ensemblGeneID$human$human_scSEG,
                                                 ctl = rownames(levinPBMCbenchmark_celSeq),
                                                 cell_type = levinPBMCbenchmark_celSeq$CellType,
                                                 # cell_type_match = T,
                                                 exprs = _,
                                                 kmeansK = rep(length(unique(levinPBMCbenchmark_celSeq$CellType)),length(unique(levinPBMCbenchmark_celSeq$batch))),
                                                 assay_name = "scMerge",
                                                 fast_svd = T))

levinPBMCbenchmark_celSeq <- runTSNE(levinPBMCbenchmark_celSeq)
plotTSNE(levinPBMCbenchmark_celSeq, colour_by = "Experiment")

levinPBMCbenchmark_celSeq <- runTSNE(levinPBMCbenchmark_celSeq, exprs_values = "scMerge")
plotTSNE(levinPBMCbenchmark_celSeq, colour_by = "Experiment")
plotTSNE(levinPBMCbenchmark_celSeq, colour_by = "CellType")
# subsetIdx <- cellTypes_pbmc10k%in%c("CD4+ T Helper2", "CD4+/CD25 T Reg", "CD4+/CD45RA+/CD25- Naive T", "CD4+/CD45RO+ Memory", "CD8+ Cytotoxic T", "CD8+/CD45RA+ Naive Cytotoxic")



levinPBMCbenchmark_celSeq$CellType2 <- as.character(droplevels(levinPBMCbenchmark_celSeq$CellType))
levinPBMCbenchmark_celSeq$CellType2[grep("monocyte", levinPBMCbenchmark_celSeq$CellType2)] <- "monocyte"
levinPBMCbenchmark_celSeq$CellType2[grep("T cell", levinPBMCbenchmark_celSeq$CellType2)] <- "T cell"




#################################################################################
######################### level 2 ###############################################
#################################################################################



set.seed(2019)



exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_celSeq, "scMerge"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_celSeq$CellType)
table(cellTypes_subset)


library(pbmcapply)

n_list <- seq(20, 400, 20)

res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = F,  n = n_list[[i]], balance = F, num_test = 100)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 50))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_celseq_level2", res_sub)


#################################################################################
######################### level 1 ###############################################
#################################################################################


exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_celSeq, "scMerge"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- (levinPBMCbenchmark_celSeq$CellType2)
table(cellTypes_subset)

set.seed(2019)



library(pbmcapply)

n_list <- seq(20, 400, 20)

res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = F,  n = n_list[[i]], balance = F, num_test = 100)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 50))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_celseq_level1", res_sub)


rm("levinPBMCbenchmark_celSeq")
gc(reset = T)



################################################
########      seqWells           ###############
################################################

levinPBMCbenchmark_seqWells <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "Seq-Well"]

levinPBMCbenchmark_seqWells$batch <- droplevels(levinPBMCbenchmark_seqWells$Experiment)

logcounts(levinPBMCbenchmark_seqWells) <- as.matrix(logcounts(levinPBMCbenchmark_seqWells))
counts(levinPBMCbenchmark_seqWells) <- as.matrix(counts(levinPBMCbenchmark_seqWells))

levinPBMCbenchmark_seqWells <- levinPBMCbenchmark_seqWells[rowSums(counts(levinPBMCbenchmark_seqWells))!=0,]
levinPBMCbenchmark_seqWells 

system.time(levinPBMCbenchmark_seqWells <- scMerge(levinPBMCbenchmark_seqWells,
                                                   # ctl = segList_ensemblGeneID$human$human_scSEG,
                                                   ctl = rownames(levinPBMCbenchmark_seqWells),
                                                   cell_type = levinPBMCbenchmark_seqWells$CellType,
                                                   # cell_type_match = T,
                                                   exprs = "logcounts",
                                                   kmeansK = rep(length(unique(levinPBMCbenchmark_seqWells$CellType)),length(unique(levinPBMCbenchmark_seqWells$batch))),
                                                   assay_name = "scMerge",
                                                   fast_svd = T))

levinPBMCbenchmark_seqWells <- runTSNE(levinPBMCbenchmark_seqWells)
plotTSNE(levinPBMCbenchmark_seqWells, colour_by = "Experiment")

levinPBMCbenchmark_seqWells <- runTSNE(levinPBMCbenchmark_seqWells, exprs_values = "scMerge")
plotTSNE(levinPBMCbenchmark_seqWells, colour_by = "Experiment")
plotTSNE(levinPBMCbenchmark_seqWells, colour_by = "CellType")
# subsetIdx <- cellTypes_pbmc10k%in%c("CD4+ T Helper2", "CD4+/CD25 T Reg", "CD4+/CD45RA+/CD25- Naive T", "CD4+/CD45RO+ Memory", "CD8+ Cytotoxic T", "CD8+/CD45RA+ Naive Cytotoxic")




levinPBMCbenchmark_seqWells$CellType2 <- as.character(droplevels(levinPBMCbenchmark_seqWells$CellType))
levinPBMCbenchmark_seqWells$CellType2[grep("monocyte", levinPBMCbenchmark_seqWells$CellType2)] <- "monocyte"
levinPBMCbenchmark_seqWells$CellType2[grep("T cell", levinPBMCbenchmark_seqWells$CellType2)] <- "T cell"


#################################################################################
######################### level 2 ###############################################
#################################################################################



exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_seqWells, "scMerge"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_seqWells$CellType)
table(cellTypes_subset)


set.seed(2019)



library(pbmcapply)

n_list <- c(seq(20, 200, 20), seq(400, 3000, 200))

res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = T,  n = n_list[[i]], balance = F, num_test = 500)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 25))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_seqWells_level2", res_sub)

#################################################################################
######################### level 1 ###############################################
#################################################################################



exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_seqWells, "scMerge"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- (levinPBMCbenchmark_seqWells$CellType2)
table(cellTypes_subset)


set.seed(2019)



library(pbmcapply)

n_list <- c(seq(20, 200, 20), seq(400, 3000, 200))

res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = T,  n = n_list[[i]], balance = F, num_test = 500)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 25))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_seqWells_level1", res_sub)

rm("levinPBMCbenchmark_seqWells")
gc(reset = T)

################################################
########      inDrops           ###############
################################################

levinPBMCbenchmark_inDrops <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "inDrops"]

levinPBMCbenchmark_inDrops$batch <- droplevels(levinPBMCbenchmark_inDrops$Experiment)

logcounts(levinPBMCbenchmark_inDrops) <- as.matrix(logcounts(levinPBMCbenchmark_inDrops))
counts(levinPBMCbenchmark_inDrops) <- as.matrix(counts(levinPBMCbenchmark_inDrops))

levinPBMCbenchmark_inDrops <- levinPBMCbenchmark_inDrops[rowSums(counts(levinPBMCbenchmark_inDrops))!=0,]
levinPBMCbenchmark_inDrops 

system.time(levinPBMCbenchmark_inDrops <- scMerge(levinPBMCbenchmark_inDrops,
                                                  # ctl = segList_ensemblGeneID$human$human_scSEG,
                                                  ctl = rownames(levinPBMCbenchmark_inDrops),
                                                  cell_type = levinPBMCbenchmark_inDrops$CellType,
                                                  # cell_type_match = T,
                                                  exprs = "logcounts",
                                                  assay_name = "scMerge",
                                                  replicate_prop = 0.25,
                                                  fast_svd = T))

levinPBMCbenchmark_inDrops <- runTSNE(levinPBMCbenchmark_inDrops)
plotTSNE(levinPBMCbenchmark_inDrops, colour_by = "Experiment")

levinPBMCbenchmark_inDrops <- runTSNE(levinPBMCbenchmark_inDrops, exprs_values = "scMerge")
plotTSNE(levinPBMCbenchmark_inDrops, colour_by = "Experiment")
plotTSNE(levinPBMCbenchmark_inDrops, colour_by = "CellType")
# subsetIdx <- cellTypes_pbmc10k%in%c("CD4+ T Helper2", "CD4+/CD25 T Reg", "CD4+/CD45RA+/CD25- Naive T", "CD4+/CD45RO+ Memory", "CD8+ Cytotoxic T", "CD8+/CD45RA+ Naive Cytotoxic")


levinPBMCbenchmark_inDrops$CellType2 <- as.character(droplevels(levinPBMCbenchmark_inDrops$CellType))
levinPBMCbenchmark_inDrops$CellType2[grep("monocyte", levinPBMCbenchmark_inDrops$CellType2)] <- "monocyte"
levinPBMCbenchmark_inDrops$CellType2[grep("T cell", levinPBMCbenchmark_inDrops$CellType2)] <- "T cell"



exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_inDrops, "scMerge"))
exprsMat_subset[exprsMat_subset<0] <- 0
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_inDrops$CellType)
table(cellTypes_subset)


system.time(res_sub <- pbmcapply::pbmclapply(1:50, function(x) {
  
  l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), 
                      cellNames = colnames(exprsMat_subset), 
                      subset_test = F,  n = 3000,  balance = F)
  table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)
}, mc.cores = 25))

rowMeans(do.call(cbind, res_sub))

#################################################################################
######################### level 2 ###############################################
#################################################################################


exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_inDrops, "scMerge"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_inDrops$CellType)
table(cellTypes_subset)


set.seed(2019)


library(pbmcapply)

n_list <- c(seq(20, 200, 20), seq(400, 5000, 200))

res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = F,  n = n_list[[i]], balance = F, num_test = 500)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 25))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_inDrops_level2", res_sub)


#################################################################################
######################### level 1 ###############################################
#################################################################################



exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_inDrops, "scMerge"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- (levinPBMCbenchmark_inDrops$CellType2)
table(cellTypes_subset)

set.seed(2019)


library(pbmcapply)

n_list <- c(seq(20, 200, 20), seq(400, 5000, 200))

res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = F,  n = n_list[[i]], balance = F, num_test = 500)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 25))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_inDrops_level1", res_sub)




rm("levinPBMCbenchmark_inDrops")
gc(reset = T)



################################################
########      dropSeqs           ###############
################################################

levinPBMCbenchmark_dropSeqs <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "Drop-seq"]

levinPBMCbenchmark_dropSeqs$batch <- droplevels(levinPBMCbenchmark_dropSeqs$Experiment)

logcounts(levinPBMCbenchmark_dropSeqs) <- as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
counts(levinPBMCbenchmark_dropSeqs) <- as.matrix(counts(levinPBMCbenchmark_dropSeqs))

levinPBMCbenchmark_dropSeqs <- levinPBMCbenchmark_dropSeqs[rowSums(counts(levinPBMCbenchmark_dropSeqs))!=0,]
levinPBMCbenchmark_dropSeqs 

system.time(levinPBMCbenchmark_dropSeqs <- scMerge(levinPBMCbenchmark_dropSeqs,
                                                   # ctl = segList_ensemblGeneID$human$human_scSEG,
                                                   ctl = rownames(levinPBMCbenchmark_dropSeqs),
                                                   cell_type = levinPBMCbenchmark_dropSeqs$CellType,
                                                   # cell_type_match = T,
                                                   exprs = "logcounts",
                                                   kmeansK = rep(length(unique(levinPBMCbenchmark_dropSeqs$CellType)),length(unique(levinPBMCbenchmark_dropSeqs$batch))),
                                                   assay_name = "scMerge",
                                                   fast_svd = T))

levinPBMCbenchmark_dropSeqs <- runTSNE(levinPBMCbenchmark_dropSeqs)
plotTSNE(levinPBMCbenchmark_dropSeqs, colour_by = "Experiment")

levinPBMCbenchmark_dropSeqs <- runTSNE(levinPBMCbenchmark_dropSeqs, exprs_values = "scMerge")
plotTSNE(levinPBMCbenchmark_dropSeqs, colour_by = "Experiment")
plotTSNE(levinPBMCbenchmark_dropSeqs, colour_by = "CellType")
# subsetIdx <- cellTypes_pbmc10k%in%c("CD4+ T Helper2", "CD4+/CD25 T Reg", "CD4+/CD45RA+/CD25- Naive T", "CD4+/CD45RO+ Memory", "CD8+ Cytotoxic T", "CD8+/CD45RA+ Naive Cytotoxic")



levinPBMCbenchmark_dropSeqs$CellType2 <- as.character(droplevels(levinPBMCbenchmark_dropSeqs$CellType))
levinPBMCbenchmark_dropSeqs$CellType2[grep("monocyte", levinPBMCbenchmark_dropSeqs$CellType2)] <- "monocyte"
levinPBMCbenchmark_dropSeqs$CellType2[grep("T cell", levinPBMCbenchmark_dropSeqs$CellType2)] <- "T cell"



#################################################################################
######################### level 2 ###############################################
#################################################################################



exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_dropSeqs, "scMerge"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_dropSeqs$CellType)
table(cellTypes_subset)


set.seed(2019)




library(pbmcapply)

n_list <- c(seq(20, 200, 20), seq(400, 5000, 200))

res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = F,  n = n_list[[i]], balance = F, num_test = 500)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 25))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_dropSeqs_level2", res_sub)




#################################################################################
######################### level 1 ###############################################
#################################################################################



exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_dropSeqs, "scMerge"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- (levinPBMCbenchmark_dropSeqs$CellType2)
table(cellTypes_subset)


set.seed(2019)




library(pbmcapply)

n_list <- c(seq(20, 200, 20), seq(400, 5000, 200))

res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = F,  n = n_list[[i]], balance = F, num_test = 500)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 25))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_dropSeqs_level1", res_sub)



################################################
########      10x Genomics       ###############
################################################

levinPBMCbenchmark_10x <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "10x Chromium (v3)"]

levinPBMCbenchmark_10x$batch <- droplevels(levinPBMCbenchmark_10x$Experiment)

logcounts(levinPBMCbenchmark_10x) <- as.matrix(logcounts(levinPBMCbenchmark_10x))
counts(levinPBMCbenchmark_10x) <- as.matrix(counts(levinPBMCbenchmark_10x))

levinPBMCbenchmark_10x <- levinPBMCbenchmark_10x[rowSums(counts(levinPBMCbenchmark_10x))!=0,]
levinPBMCbenchmark_10x 



levinPBMCbenchmark_10x$CellType2 <- as.character(droplevels(levinPBMCbenchmark_10x$CellType))
levinPBMCbenchmark_10x$CellType2[grep("monocyte", levinPBMCbenchmark_10x$CellType2)] <- "monocyte"
levinPBMCbenchmark_10x$CellType2[grep("T cell", levinPBMCbenchmark_10x$CellType2)] <- "T cell"

exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_10x, "logcounts"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_10x$CellType)
table(cellTypes_subset)


set.seed(2019)

library(pbmcapply)

n_list <- c(seq(20, 200, 20), seq(400, 3000, 200))

res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = F,  n = n_list[[i]], balance = F)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 25))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_10x_level2", res_sub)



################################################
########      Level 1       ###############
################################################

exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_10x, "logcounts"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- (levinPBMCbenchmark_10x$CellType2)
table(cellTypes_subset)


set.seed(2019)


library(pbmcapply)

n_list <- c(20, 40, 60, 80, 100, seq(200, 1600, 100), 2000)

res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = F,  n = n_list[[i]], balance = F, num_test = 500)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 25))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_10x_level1", res_sub)

rm("levinPBMCbenchmark_10x")
gc(reset = T)


################################################
########      10x Genomics (V2)      ###########
################################################

levinPBMCbenchmark_10xV2 <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "10x Chromium (v2)"]

levinPBMCbenchmark_10xV2$batch <- droplevels(levinPBMCbenchmark_10xV2$Experiment)

logcounts(levinPBMCbenchmark_10xV2) <- as.matrix(logcounts(levinPBMCbenchmark_10xV2))
counts(levinPBMCbenchmark_10xV2) <- as.matrix(counts(levinPBMCbenchmark_10xV2))

levinPBMCbenchmark_10xV2 <- levinPBMCbenchmark_10xV2[rowSums(counts(levinPBMCbenchmark_10xV2))!=0,]
levinPBMCbenchmark_10xV2 


levinPBMCbenchmark_10xV2$CellType2 <- as.character(droplevels(levinPBMCbenchmark_10xV2$CellType))
levinPBMCbenchmark_10xV2$CellType2[grep("monocyte", levinPBMCbenchmark_10xV2$CellType2)] <- "monocyte"
levinPBMCbenchmark_10xV2$CellType2[grep("T cell", levinPBMCbenchmark_10xV2$CellType2)] <- "T cell"


exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_10xV2, "logcounts"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_10xV2$CellType)
table(cellTypes_subset)


set.seed(2019)
pivot_indes <- sample(ncol(exprsMat_subset), round(ncol(exprsMat_subset)*0.5))
exprsMat_pilot <- exprsMat_subset[,pivot_indes]
dim(exprsMat_pilot)
cellTypes_pilot <- cellTypes_subset[pivot_indes]
table(cellTypes_pilot)
# 
system.time(res_sub <- mclapply(1:10, function(x) {
  
  l <- subSamplingCal(exprsMat_pilot, cellTypes_pilot, geneNames = rownames(exprsMat_pilot), cellNames = colnames(exprsMat_pilot), subset_test = F,  n = 20,  balance = F)
  table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)
}, mc.cores = 25))
gc()
# 



library(pbmcapply)


n_list <- c(20, 40, 60, 80, 100, seq(200, 1600, 100), 2000)
res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = F,  n = n_list[[i]], balance = F, num_test = 1000)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 25))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_10xV2_level2", res_sub)

################################################
########      Level 1                ###########
################################################

exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_10xV2, "logcounts"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- (levinPBMCbenchmark_10xV2$CellType2)
table(cellTypes_subset)


set.seed(2019)



library(pbmcapply)

n_list <- c(20, 40, 60, 80, 100, seq(200, 1600, 100), 2000)

res_sub <- list()
for(i in 1:length(n_list)){ 
  print(paste("n=",n_list[i]))
  system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:50, function(x) {
    tryCatch({l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), cellNames = colnames(exprsMat_subset), 
                                  subset_test = F,  n = n_list[[i]], balance = F, num_test = 500)
    table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
    error = function(e){NULL})
  }, mc.cores = 10))
  print(do.call(cbind, res_sub[[i]]))
  print(paste("Median accuracy", median(unlist(lapply(res_sub[[i]], "[[", "correct")))))
  print(paste("Median misclassified", median(unlist(lapply(res_sub[[i]], "[[", "misclassified")))))
  gc()
}


names(res_sub) <- n_list

saveRDSwithDate("sampleCalRes/levin_pbmc_10xV2_level1", res_sub)


