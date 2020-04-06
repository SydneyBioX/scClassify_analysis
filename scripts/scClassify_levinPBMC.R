
set.seed(1)
library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(CHETAH)
# library(moon)
library(pheatmap)

library(limma)
library(hopach)
library(ggraph)
library(igraph)
library(mixtools)
library(pbmcapply)
library(proxy)
library(minpack.lm)
library(sctransform)
# source("scClassify_ensembl.R")
rscript <- list.files("scClassify/")
sapply(rscript, function(x){
  source(paste("scClassify/", x, sep = ""))
})



saveRDSwithDate <- function(fileName, ...) {
  today <- gsub("-","", Sys.Date())
  fileName <- paste(paste(fileName, today, sep = "_"), ".rds", sep = "")
  saveRDS(..., file = fileName)
}



levinPBMCbenchmark <- readRDS("data/levinPBMCbenchmark.rds")


levinPBMCbenchmark_smartSeq <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "Smart-seq2"]
levinPBMCbenchmark_celSeq <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "CEL-Seq2"]
levinPBMCbenchmark_inDrops <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "inDrops"]
levinPBMCbenchmark_dropSeqs <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "Drop-seq"]
levinPBMCbenchmark_seqWells <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "Seq-Well"]
levinPBMCbenchmark_10x <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "10x Chromium (v3)"]
levinPBMCbenchmark_10xV2 <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "10x Chromium (v2)"]



levinPBMCbenchmark_smartSeq$CellType2 <- as.character(droplevels(levinPBMCbenchmark_smartSeq$CellType))
levinPBMCbenchmark_smartSeq$CellType2[grep("monocyte", levinPBMCbenchmark_smartSeq$CellType2)] <- "monocyte"
levinPBMCbenchmark_smartSeq$CellType2[grep("T cell", levinPBMCbenchmark_smartSeq$CellType2)] <- "T cell"

levinPBMCbenchmark_celSeq$CellType2 <- as.character(droplevels(levinPBMCbenchmark_celSeq$CellType))
levinPBMCbenchmark_celSeq$CellType2[grep("monocyte", levinPBMCbenchmark_celSeq$CellType2)] <- "monocyte"
levinPBMCbenchmark_celSeq$CellType2[grep("T cell", levinPBMCbenchmark_celSeq$CellType2)] <- "T cell"

levinPBMCbenchmark_10x$CellType2 <- as.character(droplevels(levinPBMCbenchmark_10x$CellType))
levinPBMCbenchmark_10x$CellType2[grep("monocyte", levinPBMCbenchmark_10x$CellType2)] <- "monocyte"
levinPBMCbenchmark_10x$CellType2[grep("T cell", levinPBMCbenchmark_10x$CellType2)] <- "T cell"

levinPBMCbenchmark_10xV2$CellType2 <- as.character(droplevels(levinPBMCbenchmark_10xV2$CellType))
levinPBMCbenchmark_10xV2$CellType2[grep("monocyte", levinPBMCbenchmark_10xV2$CellType2)] <- "monocyte"
levinPBMCbenchmark_10xV2$CellType2[grep("T cell", levinPBMCbenchmark_10xV2$CellType2)] <- "T cell"

levinPBMCbenchmark_seqWells$CellType2 <- as.character(droplevels(levinPBMCbenchmark_seqWells$CellType))
levinPBMCbenchmark_seqWells$CellType2[grep("monocyte", levinPBMCbenchmark_seqWells$CellType2)] <- "monocyte"
levinPBMCbenchmark_seqWells$CellType2[grep("T cell", levinPBMCbenchmark_seqWells$CellType2)] <- "T cell"

levinPBMCbenchmark_inDrops$CellType2 <- as.character(droplevels(levinPBMCbenchmark_inDrops$CellType))
levinPBMCbenchmark_inDrops$CellType2[grep("monocyte", levinPBMCbenchmark_inDrops$CellType2)] <- "monocyte"
levinPBMCbenchmark_inDrops$CellType2[grep("T cell", levinPBMCbenchmark_inDrops$CellType2)] <- "T cell"

levinPBMCbenchmark_dropSeqs$CellType2 <- as.character(droplevels(levinPBMCbenchmark_dropSeqs$CellType))
levinPBMCbenchmark_dropSeqs$CellType2[grep("monocyte", levinPBMCbenchmark_dropSeqs$CellType2)] <- "monocyte"
levinPBMCbenchmark_dropSeqs$CellType2[grep("T cell", levinPBMCbenchmark_dropSeqs$CellType2)] <- "T cell"
# logcounts(levinPBMCbenchmark_smartSeq)



set.seed(1)
system.time(trainInDropsRes_prob70 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                 cellTypes_train = levinPBMCbenchmark_inDrops$CellType2,
                                                 exprsMat_test = list(
                                                   smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                   celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                   tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                   tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                   seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                   inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                   dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                 ),
                                                 cellTypes_test = list(
                                                   smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                                   celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                                   tenX = levinPBMCbenchmark_10x$CellType2,
                                                   tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                                   seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                                   inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                                   dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                                 similarity = c("pearson"),
                                                 selectFeatures = c("limma"),
                                                 algorithm = c("WKNN"),
                                                 parallel = F,
                                                 ncores = 1,
                                                 prob_threshold = 0.7))
do.call(rbind, lapply(trainInDropsRes_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

set.seed(1)
system.time(trainSeqWellsRes_prob70 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                  cellTypes_train = levinPBMCbenchmark_seqWells$CellType2,
                                                  exprsMat_test = list(
                                                    smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                    celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                    tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                    tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                    seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                    inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                    dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                  ),
                                                  cellTypes_test = list(
                                                    smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                                    celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                                    tenX = levinPBMCbenchmark_10x$CellType2,
                                                    tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                                    seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                                    inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                                    dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                                  similarity = c("pearson"),
                                                  selectFeatures = c("limma"),
                                                  algorithm = c("WKNN"),
                                                  parallel = F,
                                                  ncores = 1,
                                                  prob_threshold = 0.7))
do.call(rbind, lapply(trainSeqWellsRes_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

set.seed(1)
system.time(trainDropSeqsRes_prob70 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs)),
                                                  cellTypes_train = levinPBMCbenchmark_dropSeqs$CellType2,
                                                  exprsMat_test = list(
                                                    smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                    celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                    tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                    tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                    seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                    inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                    dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                  ),
                                                  cellTypes_test = list(
                                                    smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                                    celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                                    tenX = levinPBMCbenchmark_10x$CellType2,
                                                    tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                                    seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                                    inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                                    dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                                  similarity = c("pearson"),
                                                  selectFeatures = c("limma"),
                                                  algorithm = c("WKNN"),
                                                  parallel = F,
                                                  ncores = 1,
                                                  prob_threshold = 0.7))
do.call(rbind, lapply(trainDropSeqsRes_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

set.seed(1)
system.time(train10xRes_prob70 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                             cellTypes_train = levinPBMCbenchmark_10x$CellType2,
                                             exprsMat_test = list(
                                               smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                               celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                               tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                               tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                               seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                               inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                               dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                             ),
                                             cellTypes_test = list(
                                               smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                               celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                               tenX = levinPBMCbenchmark_10x$CellType2,
                                               tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                               seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                               inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                               dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                             similarity = c("pearson"),
                                             selectFeatures = c("limma"),
                                             algorithm = c("WKNN"),
                                             parallel = F,
                                             ncores = 1,
                                             prob_threshold = 0.7))
do.call(rbind, lapply(train10xRes_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

set.seed(1)
system.time(train10xV2Res_prob70 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                               cellTypes_train = levinPBMCbenchmark_10xV2$CellType2,
                                               exprsMat_test = list(
                                                 smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                 celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                 tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                 tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                 seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                 inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                 dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                               ),
                                               cellTypes_test = list(
                                                 smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                                 celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                                 tenX = levinPBMCbenchmark_10x$CellType2,
                                                 tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                                 seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                                 inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                                 dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                               similarity = c("pearson"),
                                               selectFeatures = c("limma"),
                                               algorithm = c("WKNN"),
                                               parallel = F,
                                               ncores = 1,
                                               prob_threshold = 0.7))
do.call(rbind, lapply(train10xV2Res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))



set.seed(1)
system.time(trainSmartSeqRes_prob70 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                  cellTypes_train = levinPBMCbenchmark_smartSeq$CellType2,
                                                  exprsMat_test = list(
                                                    smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                    celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                    tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                    tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                    seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                    inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                    dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                  ),
                                                  cellTypes_test = list(
                                                    smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                                    celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                                    tenX = levinPBMCbenchmark_10x$CellType2,
                                                    tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                                    seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                                    inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                                    dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                                  similarity = c("pearson"),
                                                  selectFeatures = c("limma"),
                                                  algorithm = c("WKNN"),
                                                  parallel = F,
                                                  ncores = 5,
                                                  prob_threshold = 0.7))

do.call(rbind, lapply(trainSmartSeqRes_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))


set.seed(1)
system.time(trainCelSeqRes_prob70 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                cellTypes_train = levinPBMCbenchmark_celSeq$CellType2,
                                                exprsMat_test = list(
                                                  smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                  celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                  tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                  tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                  seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                  inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                  dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                ),
                                                cellTypes_test = list(
                                                  smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                                  celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                                  tenX = levinPBMCbenchmark_10x$CellType2,
                                                  tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                                  seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                                  inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                                  dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                                similarity = c("pearson"),
                                                selectFeatures = c("limma"),
                                                algorithm = c("WKNN"),
                                                parallel = F,
                                                ncores = 1,
                                                prob_threshold = 0.7))

do.call(rbind, lapply(trainCelSeqRes_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

scClassifyRes_prob70 <- rbind(
  do.call(rbind, lapply(train10xRes_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(train10xV2Res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainCelSeqRes_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainDropSeqsRes_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainInDropsRes_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainSeqWellsRes_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainSmartSeqRes_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))
)

rownames(scClassifyRes_prob70) <- paste(c(rep("tenX", 7),
                                          rep("tenXV2", 7),
                                          rep("celSeq", 7), 
                                          rep("dropSeqs", 7),
                                          rep("inDrops", 7),
                                          rep("seqWells", 7),
                                          rep("smartseq", 7)), rownames(scClassifyRes_prob70), sep = "_")

scClassifyRes_prob70 <- scClassifyRes_prob70[unlist(lapply(strsplit(rownames(scClassifyRes_prob70), "_"), function(x) x[1] != x[2])),]
saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/results/levinPBMCResults/levinPBMC_scClassifyRes_prob70", scClassifyRes_prob70)




#################################################################################################################################################################################################################################
################################################################################ Ensemble ###################################################################################################################
############################################################################################################################################################################################################


set.seed(1)
system.time(trainSmartSeqRes_prob70_ensemble <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                           cellTypes_train = levinPBMCbenchmark_smartSeq$CellType2,
                                                           exprsMat_test = list(
                                                             smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                             celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                             tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                             tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                             seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                             inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                             dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                           ),
                                                           cellTypes_test = list(
                                                             smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                                             celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                                             tenX = levinPBMCbenchmark_10x$CellType2,
                                                             tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                                             seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                                             inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                                             dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                                           similarity = c("pearson",  "spearman",
                                                                          "cosine", "jaccard", "kendall",
                                                                          "weighted_rank"),
                                                           selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                           algorithm = c("WKNN"),
                                                           parallel = T,
                                                           ncores = 15,
                                                           prob_threshold = 0.7))

saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinSmartSeqRes_prob70_ensemble", trainSmartSeqRes_prob70_ensemble)



set.seed(1)
system.time(trainCelSeqRes_prob70_ensemble <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                         cellTypes_train = levinPBMCbenchmark_celSeq$CellType2,
                                                             exprsMat_test = list(
                                                               smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                               celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                               tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                               tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                               seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                               inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                               dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                             ),
                                                             cellTypes_test = list(
                                                               smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                                               celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                                               tenX = levinPBMCbenchmark_10x$CellType2,
                                                               tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                                               seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                                               inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                                               dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                                             similarity = c("pearson",  "spearman",
                                                                            "cosine", "jaccard", "kendall",
                                                                            "weighted_rank"),
                                                             selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                             algorithm = c("WKNN", "KNN"),
                                                             parallel = T,
                                                             ncores = 5,
                                                             prob_threshold = 0.7))

saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinCelSeqRes_prob70_ensemble", trainCelSeqRes_prob70_ensemble)


set.seed(1)
system.time(trainInDropsRes_prob70_ensemble <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                 cellTypes_train = levinPBMCbenchmark_inDrops$CellType2,
                                                 exprsMat_test = list(
                                                   smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                   celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                   tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                   tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                   seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                   inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                   dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                 ),
                                                 cellTypes_test = list(
                                                   smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                                   celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                                   tenX = levinPBMCbenchmark_10x$CellType2,
                                                   tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                                   seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                                   inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                                   dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                                 similarity = c("pearson",  "spearman",
                                                                "cosine", "jaccard", "kendall",
                                                                "weighted_rank"),
                                                 selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                 algorithm = c("WKNN"),
                                                 parallel = T,
                                                 ncores = 15,
                                                 prob_threshold = 0.7))
saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinInDropsRes_prob70_ensemble", trainInDropsRes_prob70_ensemble)

set.seed(1)
system.time(trainSeqWellsRes_prob70_ensemble <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                  cellTypes_train = levinPBMCbenchmark_seqWells$CellType2,
                                                  exprsMat_test = list(
                                                    smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                    celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                    tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                    tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                    seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                    inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                    dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                  ),
                                                  cellTypes_test = list(
                                                    smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                                    celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                                    tenX = levinPBMCbenchmark_10x$CellType2,
                                                    tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                                    seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                                    inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                                    dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                                  similarity = c("pearson",  "spearman",
                                                                 "cosine", "jaccard", "kendall",
                                                                 "weighted_rank"),
                                                  selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                  algorithm = c("WKNN"),
                                                  parallel = T,
                                                  ncores = 15,
                                                  prob_threshold = 0.7))
saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinSeqWellsRes_prob70_ensemble", trainSeqWellsRes_prob70_ensemble)



set.seed(1)
system.time(trainDropSeqsRes_prob70_ensemble <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs)),
                                                  cellTypes_train = levinPBMCbenchmark_dropSeqs$CellType2,
                                                  exprsMat_test = list(
                                                    smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                    celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                    tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                    tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                    seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                    inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                    dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                  ),
                                                  cellTypes_test = list(
                                                    smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                                    celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                                    tenX = levinPBMCbenchmark_10x$CellType2,
                                                    tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                                    seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                                    inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                                    dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                                  similarity = c("pearson",  "spearman",
                                                                 "cosine", "jaccard", "kendall",
                                                                 "weighted_rank"),
                                                  selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                  algorithm = c("WKNN"),
                                                  parallel = T,
                                                  ncores = 15,
                                                  prob_threshold = 0.7))
saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinDropSeqsRes_prob70_ensemble", trainDropSeqsRes_prob70_ensemble)


set.seed(1)
system.time(train10xRes_prob70_ensemble <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                             cellTypes_train = levinPBMCbenchmark_10x$CellType2,
                                             exprsMat_test = list(
                                               smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                               celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                               tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                               tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                               seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                               inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                               dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                             ),
                                             cellTypes_test = list(
                                               smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                               celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                               tenX = levinPBMCbenchmark_10x$CellType2,
                                               tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                               seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                               inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                               dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                             similarity = c("pearson",  "spearman",
                                                            "cosine", "jaccard", "kendall",
                                                            "weighted_rank"),
                                             selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                             algorithm = c("WKNN", "KNN"),
                                             parallel = T,
                                             ncores = 15,
                                             prob_threshold = 0.7))
saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinTenXRes_prob70_ensemble", train10xRes_prob70_ensemble)


set.seed(1)
system.time(train10xV2Res_prob70_ensemble <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                               cellTypes_train = levinPBMCbenchmark_10xV2$CellType2,
                                               exprsMat_test = list(
                                                 smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                 celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                 tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                 tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                 seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                 inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                 dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                               ),
                                               cellTypes_test = list(
                                                 smartseq = levinPBMCbenchmark_smartSeq$CellType2,
                                                 celSeq = levinPBMCbenchmark_celSeq$CellType2,
                                                 tenX = levinPBMCbenchmark_10x$CellType2,
                                                 tenXV2 = levinPBMCbenchmark_10xV2$CellType2,
                                                 seqWells = levinPBMCbenchmark_seqWells$CellType2,
                                                 inDrops = levinPBMCbenchmark_inDrops$CellType2,
                                                 dropSeqs = levinPBMCbenchmark_dropSeqs$CellType2),
                                               similarity = c("pearson",  "spearman",
                                                              "cosine", "jaccard", "kendall",
                                                              "weighted_rank"),
                                               selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                               algorithm = c("WKNN"),
                                               parallel = T,
                                               ncores = 5,
                                               prob_threshold = 0.7))
saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinTenXV2Res_prob70_ensemble", train10xV2Res_prob70_ensemble)


###################################################################################################################################################
###################################################################################################################################################
###############################################   Level 2      ####################################################################################
###################################################################################################################################################

set.seed(1)
system.time(trainInDropsRes_prob60_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                 cellTypes_train = levinPBMCbenchmark_inDrops$CellType,
                                                 exprsMat_test = list(
                                                   smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                   celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                   tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                   tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                   seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                   inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                   dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                 ),
                                                 cellTypes_test = list(
                                                   smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                   celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                   tenX = levinPBMCbenchmark_10x$CellType,
                                                   tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                   seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                   inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                   dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                 similarity = c("pearson"),
                                                 selectFeatures = c("limma"),
                                                 algorithm = c("WKNN"),
                                                 parallel = F,
                                                 ncores = 1,
                                                 prob_threshold = 0.6))
do.call(rbind, lapply(trainInDropsRes_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

set.seed(1)
system.time(trainSeqWellsRes_prob60_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                  cellTypes_train = levinPBMCbenchmark_seqWells$CellType,
                                                  exprsMat_test = list(
                                                    smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                    celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                    tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                    tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                    seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                    inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                    dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                  ),
                                                  cellTypes_test = list(
                                                    smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                    celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                    tenX = levinPBMCbenchmark_10x$CellType,
                                                    tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                    seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                    inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                    dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                  similarity = c("pearson"),
                                                  selectFeatures = c("limma"),
                                                  algorithm = c("WKNN"),
                                                  parallel = F,
                                                  ncores = 1,
                                                  prob_threshold = 0.6))
do.call(rbind, lapply(trainSeqWellsRes_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

set.seed(1)
system.time(trainDropSeqsRes_prob60_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs)),
                                                  cellTypes_train = levinPBMCbenchmark_dropSeqs$CellType,
                                                  exprsMat_test = list(
                                                    smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                    celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                    tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                    tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                    seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                    inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                    dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                  ),
                                                  cellTypes_test = list(
                                                    smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                    celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                    tenX = levinPBMCbenchmark_10x$CellType,
                                                    tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                    seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                    inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                    dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                  similarity = c("pearson"),
                                                  selectFeatures = c("limma"),
                                                  algorithm = c("WKNN"),
                                                  parallel = F,
                                                  ncores = 1,
                                                  prob_threshold = 0.6))
do.call(rbind, lapply(trainDropSeqsRes_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

set.seed(1)
system.time(train10xRes_prob60_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                             cellTypes_train = levinPBMCbenchmark_10x$CellType,
                                             exprsMat_test = list(
                                               smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                               celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                               tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                               tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                               seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                               inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                               dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                             ),
                                             cellTypes_test = list(
                                               smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                               celSeq = levinPBMCbenchmark_celSeq$CellType,
                                               tenX = levinPBMCbenchmark_10x$CellType,
                                               tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                               seqWells = levinPBMCbenchmark_seqWells$CellType,
                                               inDrops = levinPBMCbenchmark_inDrops$CellType,
                                               dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                             similarity = c("pearson"),
                                             selectFeatures = c("limma"),
                                             algorithm = c("WKNN"),
                                             parallel = F,
                                             ncores = 1,
                                             prob_threshold = 0.6))
do.call(rbind, lapply(train10xRes_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

set.seed(1)
system.time(train10xV2Res_prob60_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                               cellTypes_train = levinPBMCbenchmark_10xV2$CellType,
                                               exprsMat_test = list(
                                                 smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                 celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                 tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                 tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                 seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                 inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                 dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                               ),
                                               cellTypes_test = list(
                                                 smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                 celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                 tenX = levinPBMCbenchmark_10x$CellType,
                                                 tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                 seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                 inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                 dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                               similarity = c("pearson"),
                                               selectFeatures = c("limma"),
                                               algorithm = c("WKNN"),
                                               parallel = F,
                                               ncores = 1,
                                               prob_threshold = 0.6))
do.call(rbind, lapply(train10xV2Res_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))



set.seed(1)
system.time(trainSmartSeqRes_prob60_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                  cellTypes_train = levinPBMCbenchmark_smartSeq$CellType,
                                                  exprsMat_test = list(
                                                    smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                    celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                    tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                    tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                    seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                    inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                    dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                  ),
                                                  cellTypes_test = list(
                                                    smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                    celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                    tenX = levinPBMCbenchmark_10x$CellType,
                                                    tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                    seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                    inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                    dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                  similarity = c("pearson"),
                                                  selectFeatures = c("limma"),
                                                  algorithm = c("WKNN"),
                                                  parallel = F,
                                                  ncores = 5,
                                                  prob_threshold = 0.6))

do.call(rbind, lapply(trainSmartSeqRes_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))


set.seed(1)
system.time(trainCelSeqRes_prob60_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                cellTypes_train = levinPBMCbenchmark_celSeq$CellType,
                                                exprsMat_test = list(
                                                  smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                  celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                  tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                  tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                  seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                  inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                  dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                ),
                                                cellTypes_test = list(
                                                  smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                  celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                  tenX = levinPBMCbenchmark_10x$CellType,
                                                  tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                  seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                  inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                  dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                similarity = c("pearson"),
                                                selectFeatures = c("limma"),
                                                algorithm = c("WKNN"),
                                                parallel = F,
                                                ncores = 1,
                                                prob_threshold = 0.6))

do.call(rbind, lapply(trainCelSeqRes_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

scClassifyRes_prob60_level2 <- rbind(
  do.call(rbind, lapply(train10xRes_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(train10xV2Res_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainCelSeqRes_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainDropSeqsRes_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainInDropsRes_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainSeqWellsRes_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainSmartSeqRes_prob60_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))
)

rownames(scClassifyRes_prob60_level2) <- paste(c(rep("tenX", 7),
                                          rep("tenXV2", 7),
                                          rep("celSeq", 7), 
                                          rep("dropSeqs", 7),
                                          rep("inDrops", 7),
                                          rep("seqWells", 7),
                                          rep("smartseq", 7)), rownames(scClassifyRes_prob60_level2), sep = "_")

scClassifyRes_prob60_level2 <- scClassifyRes_prob60_level2[unlist(lapply(strsplit(rownames(scClassifyRes_prob60_level2), "_"), function(x) x[1] != x[2])),]
saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/results/levinPBMCResults/levinPBMC_scClassifyRes_prob60_level2", scClassifyRes_prob60_level2)




set.seed(1)
system.time(trainInDropsRes_prob70_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                        cellTypes_train = levinPBMCbenchmark_inDrops$CellType,
                                                        exprsMat_test = list(
                                                          smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                          celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                          tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                          tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                          seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                          inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                          dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                        ),
                                                        cellTypes_test = list(
                                                          smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                          celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                          tenX = levinPBMCbenchmark_10x$CellType,
                                                          tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                          seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                          inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                          dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                        similarity = c("pearson"),
                                                        selectFeatures = c("limma"),
                                                        algorithm = c("WKNN"),
                                                        parallel = F,
                                                        ncores = 1,
                                                        prob_threshold = 0.7))
do.call(rbind, lapply(trainInDropsRes_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

set.seed(1)
system.time(trainSeqWellsRes_prob70_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                         cellTypes_train = levinPBMCbenchmark_seqWells$CellType,
                                                         exprsMat_test = list(
                                                           smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                           celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                           tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                           tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                           seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                           inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                           dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                         ),
                                                         cellTypes_test = list(
                                                           smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                           celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                           tenX = levinPBMCbenchmark_10x$CellType,
                                                           tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                           seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                           inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                           dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                         similarity = c("pearson"),
                                                         selectFeatures = c("limma"),
                                                         algorithm = c("WKNN"),
                                                         parallel = F,
                                                         ncores = 1,
                                                         prob_threshold = 0.7))
do.call(rbind, lapply(trainSeqWellsRes_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

set.seed(1)
system.time(trainDropSeqsRes_prob70_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs)),
                                                         cellTypes_train = levinPBMCbenchmark_dropSeqs$CellType,
                                                         exprsMat_test = list(
                                                           smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                           celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                           tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                           tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                           seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                           inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                           dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                         ),
                                                         cellTypes_test = list(
                                                           smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                           celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                           tenX = levinPBMCbenchmark_10x$CellType,
                                                           tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                           seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                           inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                           dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                         similarity = c("pearson"),
                                                         selectFeatures = c("limma"),
                                                         algorithm = c("WKNN"),
                                                         parallel = F,
                                                         ncores = 1,
                                                         prob_threshold = 0.7))
do.call(rbind, lapply(trainDropSeqsRes_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

set.seed(1)
system.time(train10xRes_prob70_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                    cellTypes_train = levinPBMCbenchmark_10x$CellType,
                                                    exprsMat_test = list(
                                                      smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                      celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                      tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                      tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                      seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                      inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                      dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                    ),
                                                    cellTypes_test = list(
                                                      smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                      celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                      tenX = levinPBMCbenchmark_10x$CellType,
                                                      tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                      seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                      inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                      dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                    similarity = c("pearson"),
                                                    selectFeatures = c("limma"),
                                                    algorithm = c("WKNN"),
                                                    parallel = F,
                                                    ncores = 1,
                                                    prob_threshold = 0.7))
do.call(rbind, lapply(train10xRes_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

set.seed(1)
system.time(train10xV2Res_prob70_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                      cellTypes_train = levinPBMCbenchmark_10xV2$CellType,
                                                      exprsMat_test = list(
                                                        smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                        celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                        tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                        tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                        seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                        inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                        dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                      ),
                                                      cellTypes_test = list(
                                                        smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                        celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                        tenX = levinPBMCbenchmark_10x$CellType,
                                                        tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                        seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                        inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                        dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                      similarity = c("pearson"),
                                                      selectFeatures = c("limma"),
                                                      algorithm = c("WKNN"),
                                                      parallel = F,
                                                      ncores = 1,
                                                      prob_threshold = 0.7))
do.call(rbind, lapply(train10xV2Res_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))



set.seed(1)
system.time(trainSmartSeqRes_prob70_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                         cellTypes_train = levinPBMCbenchmark_smartSeq$CellType,
                                                         exprsMat_test = list(
                                                           smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                           celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                           tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                           tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                           seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                           inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                           dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                         ),
                                                         cellTypes_test = list(
                                                           smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                           celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                           tenX = levinPBMCbenchmark_10x$CellType,
                                                           tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                           seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                           inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                           dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                         similarity = c("pearson"),
                                                         selectFeatures = c("limma"),
                                                         algorithm = c("WKNN"),
                                                         parallel = F,
                                                         ncores = 5,
                                                         prob_threshold = 0.7))

do.call(rbind, lapply(trainSmartSeqRes_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))


set.seed(1)
system.time(trainCelSeqRes_prob70_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                       cellTypes_train = levinPBMCbenchmark_celSeq$CellType,
                                                       exprsMat_test = list(
                                                         smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                         celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                         tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                         tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                         seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                         inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                         dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                       ),
                                                       cellTypes_test = list(
                                                         smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                         celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                         tenX = levinPBMCbenchmark_10x$CellType,
                                                         tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                         seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                         inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                         dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                       similarity = c("pearson"),
                                                       selectFeatures = c("limma"),
                                                       algorithm = c("WKNN"),
                                                       parallel = F,
                                                       ncores = 1,
                                                       prob_threshold = 0.7))

do.call(rbind, lapply(trainCelSeqRes_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

scClassifyRes_prob70_level2 <- rbind(
  do.call(rbind, lapply(train10xRes_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(train10xV2Res_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainCelSeqRes_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainDropSeqsRes_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainInDropsRes_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainSeqWellsRes_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(trainSmartSeqRes_prob70_level2$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))
)

rownames(scClassifyRes_prob70_level2) <- paste(c(rep("tenX", 7),
                                                 rep("tenXV2", 7),
                                                 rep("celSeq", 7), 
                                                 rep("dropSeqs", 7),
                                                 rep("inDrops", 7),
                                                 rep("seqWells", 7),
                                                 rep("smartseq", 7)), rownames(scClassifyRes_prob70_level2), sep = "_")

scClassifyRes_prob70_level2 <- scClassifyRes_prob70_level2[unlist(lapply(strsplit(rownames(scClassifyRes_prob70_level2), "_"), function(x) x[1] != x[2])),]
saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/results/levinPBMCResults/levinPBMC_scClassifyRes_prob70_level2", scClassifyRes_prob70_level2)








set.seed(1)
system.time(trainSmartSeqRes_prob70_ensemble_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                           cellTypes_train = levinPBMCbenchmark_smartSeq$CellType,
                                                           exprsMat_test = list(
                                                             smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                             celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                             tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                             tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                             seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                             inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                             dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                           ),
                                                           cellTypes_test = list(
                                                             smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                             celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                             tenX = levinPBMCbenchmark_10x$CellType,
                                                             tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                             seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                             inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                             dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                           similarity = c("pearson",  "spearman",
                                                                          "cosine", "jaccard", "kendall",
                                                                          "weighted_rank"),
                                                           selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                           algorithm = c("WKNN"),
                                                           parallel = T,
                                                           ncores = 15,
                                                           prob_threshold = 0.7))

saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinSmartSeqRes_prob70_ensemble_level2", trainSmartSeqRes_prob70_ensemble_level2)



set.seed(1)
system.time(trainCelSeqRes_prob70_ensemble_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                         cellTypes_train = levinPBMCbenchmark_celSeq$CellType,
                                                         exprsMat_test = list(
                                                           smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                           celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                           tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                           tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                           seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                           inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                           dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                         ),
                                                         cellTypes_test = list(
                                                           smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                           celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                           tenX = levinPBMCbenchmark_10x$CellType,
                                                           tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                           seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                           inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                           dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                         similarity = c("pearson",  "spearman",
                                                                        "cosine", "jaccard", "kendall",
                                                                        "weighted_rank"),
                                                         selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                         algorithm = c("WKNN", "KNN"),
                                                         parallel = T,
                                                         ncores = 5,
                                                         prob_threshold = 0.7))

saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinCelSeqRes_prob70_ensemble_level2", trainCelSeqRes_prob70_ensemble_level2)


set.seed(1)
system.time(trainInDropsRes_prob70_ensemble_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                          cellTypes_train = levinPBMCbenchmark_inDrops$CellType,
                                                          exprsMat_test = list(
                                                            smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                            celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                            tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                            tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                            seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                            inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                            dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                          ),
                                                          cellTypes_test = list(
                                                            smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                            celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                            tenX = levinPBMCbenchmark_10x$CellType,
                                                            tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                            seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                            inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                            dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                          similarity = c("pearson",  "spearman",
                                                                         "cosine", "jaccard", "kendall",
                                                                         "weighted_rank"),
                                                          selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                          algorithm = c("WKNN"),
                                                          parallel = T,
                                                          ncores = 15,
                                                          prob_threshold = 0.7))
saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinInDropsRes_prob70_ensemble_level2", trainInDropsRes_prob70_ensemble_level2)

set.seed(1)
system.time(trainSeqWellsRes_prob70_ensemble_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                           cellTypes_train = levinPBMCbenchmark_seqWells$CellType,
                                                           exprsMat_test = list(
                                                             smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                             celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                             tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                             tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                             seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                             inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                             dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                           ),
                                                           cellTypes_test = list(
                                                             smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                             celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                             tenX = levinPBMCbenchmark_10x$CellType,
                                                             tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                             seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                             inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                             dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                           similarity = c("pearson",  "spearman",
                                                                          "cosine", "jaccard", "kendall",
                                                                          "weighted_rank"),
                                                           selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                           algorithm = c("WKNN"),
                                                           parallel = T,
                                                           ncores = 15,
                                                           prob_threshold = 0.7))
saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinSeqWellsRes_prob70_ensemble_level2", trainSeqWellsRes_prob70_ensemble_level2)



set.seed(1)
system.time(trainDropSeqsRes_prob70_ensemble_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs)),
                                                           cellTypes_train = levinPBMCbenchmark_dropSeqs$CellType,
                                                           exprsMat_test = list(
                                                             smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                             celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                             tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                             tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                             seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                             inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                             dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                           ),
                                                           cellTypes_test = list(
                                                             smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                             celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                             tenX = levinPBMCbenchmark_10x$CellType,
                                                             tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                             seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                             inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                             dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                           similarity = c("pearson",  "spearman",
                                                                          "cosine", "jaccard", "kendall",
                                                                          "weighted_rank"),
                                                           selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                           algorithm = c("WKNN"),
                                                           parallel = T,
                                                           ncores = 15,
                                                           prob_threshold = 0.7))
saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinDropSeqsRes_prob70_ensemble_level2", trainDropSeqsRes_prob70_ensemble_level2)


set.seed(1)
system.time(train10xRes_prob70_ensemble_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                      cellTypes_train = levinPBMCbenchmark_10x$CellType,
                                                      exprsMat_test = list(
                                                        smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                        celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                        tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                        tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                        seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                        inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                        dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                      ),
                                                      cellTypes_test = list(
                                                        smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                        celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                        tenX = levinPBMCbenchmark_10x$CellType,
                                                        tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                        seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                        inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                        dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                      similarity = c("pearson",  "spearman",
                                                                     "cosine", "jaccard", "kendall",
                                                                     "weighted_rank"),
                                                      selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                      algorithm = c("WKNN", "KNN"),
                                                      parallel = T,
                                                      ncores = 15,
                                                      prob_threshold = 0.7))
saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinTenXRes_prob70_ensemble_level2", train10xRes_prob70_ensemble_level2)


set.seed(1)
system.time(train10xV2Res_prob70_ensemble_level2 <- scClassify(exprsMat_train = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                        cellTypes_train = levinPBMCbenchmark_10xV2$CellType,
                                                        exprsMat_test = list(
                                                          smartseq = as.matrix(logcounts(levinPBMCbenchmark_smartSeq)),
                                                          celSeq = as.matrix(logcounts(levinPBMCbenchmark_celSeq)),
                                                          tenX = as.matrix(logcounts(levinPBMCbenchmark_10x)),
                                                          tenXV2 = as.matrix(logcounts(levinPBMCbenchmark_10xV2)),
                                                          seqWells = as.matrix(logcounts(levinPBMCbenchmark_seqWells)),
                                                          inDrops = as.matrix(logcounts(levinPBMCbenchmark_inDrops)),
                                                          dropSeqs = as.matrix(logcounts(levinPBMCbenchmark_dropSeqs))
                                                        ),
                                                        cellTypes_test = list(
                                                          smartseq = levinPBMCbenchmark_smartSeq$CellType,
                                                          celSeq = levinPBMCbenchmark_celSeq$CellType,
                                                          tenX = levinPBMCbenchmark_10x$CellType,
                                                          tenXV2 = levinPBMCbenchmark_10xV2$CellType,
                                                          seqWells = levinPBMCbenchmark_seqWells$CellType,
                                                          inDrops = levinPBMCbenchmark_inDrops$CellType,
                                                          dropSeqs = levinPBMCbenchmark_dropSeqs$CellType),
                                                        similarity = c("pearson",  "spearman",
                                                                       "cosine", "jaccard", "kendall",
                                                                       "weighted_rank"),
                                                        selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                        algorithm = c("WKNN"),
                                                        parallel = T,
                                                        ncores = 5,
                                                        prob_threshold = 0.7))
saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/levinTenXV2Res_prob70_ensemble_level2", train10xV2Res_prob70_ensemble_level2)




##############################################################################
#######################      Joint classification weights       ##############
##############################################################################



library(scMerge)
data("segList_ensemblGeneID")
source("sampleCal_lfunctions.R")



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


exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_smartSeq, "logcounts"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_smartSeq$CellType)
table(cellTypes_subset)




system.time(res_sub <- mclapply(1:50, function(x) {
  
  l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), 
                      cellNames = colnames(exprsMat_subset), subset_test = F,  n = round(ncol(exprsMat_subset)/5)*4,  balance = F)
  table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)
}, mc.cores = 10))
gc()


res_smartseq <- do.call(cbind, res_sub)

################################################
########      Cel-seq2         ###############
################################################


levinPBMCbenchmark_celSeq <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "CEL-Seq2"]

levinPBMCbenchmark_celSeq$batch <- droplevels(levinPBMCbenchmark_celSeq$Experiment)

logcounts(levinPBMCbenchmark_celSeq) <- as.matrix(logcounts(levinPBMCbenchmark_celSeq))
counts(levinPBMCbenchmark_celSeq) <- as.matrix(counts(levinPBMCbenchmark_celSeq))

levinPBMCbenchmark_celSeq <- levinPBMCbenchmark_celSeq[rowSums(counts(levinPBMCbenchmark_celSeq))!=0,]
levinPBMCbenchmark_celSeq 

exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_celSeq, "logcounts"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_celSeq$CellType)
table(cellTypes_subset)




system.time(res_sub <- mclapply(1:50, function(x) {
  
  l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), 
                      cellNames = colnames(exprsMat_subset), subset_test = F,  n = round(ncol(exprsMat_subset)/5)*4,  balance = F)
  table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)
}, mc.cores = 10))
gc()


res_celseq <- do.call(cbind, res_sub)


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
                                                  kmeansK = rep(length(unique(levinPBMCbenchmark_inDrops$CellType)),length(unique(levinPBMCbenchmark_inDrops$batch))),
                                                  assay_name = "scMerge",
                                                  fast_svd = T))

levinPBMCbenchmark_inDrops <- runTSNE(levinPBMCbenchmark_inDrops)
plotTSNE(levinPBMCbenchmark_inDrops, colour_by = "Experiment")

levinPBMCbenchmark_inDrops <- runTSNE(levinPBMCbenchmark_inDrops, exprs_values = "scMerge")
plotTSNE(levinPBMCbenchmark_inDrops, colour_by = "Experiment")
plotTSNE(levinPBMCbenchmark_inDrops, colour_by = "CellType")
# subsetIdx <- cellTypes_pbmc10k%in%c("CD4+ T Helper2", "CD4+/CD25 T Reg", "CD4+/CD45RA+/CD25- Naive T", "CD4+/CD45RO+ Memory", "CD8+ Cytotoxic T", "CD8+/CD45RA+ Naive Cytotoxic")

exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_inDrops, "scMerge"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_inDrops$CellType)
table(cellTypes_subset)


set.seed(2019)



system.time(res_sub <- pbmcapply::pbmclapply(1:50, function(x) {
  
  l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), 
                      cellNames = colnames(exprsMat_subset), subset_test = F,  n = round(ncol(exprsMat_subset)/5)*4,  balance = F)
  table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)
}, mc.cores = 10))
gc()


res_indrops <- do.call(cbind, res_sub)


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

exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_dropSeqs, "scMerge"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_dropSeqs$CellType)
table(cellTypes_subset)


system.time(res_sub <- pbmcapply::pbmclapply(1:50, function(x) {
  
  l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), 
                      cellNames = colnames(exprsMat_subset), subset_test = F,  n = round(ncol(exprsMat_subset)/5)*4,  balance = F)
  table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)
}, mc.cores = 10))
gc()


res_dropseq <- do.call(cbind, res_sub)





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

exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_seqWells, "scMerge"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_seqWells$CellType)
table(cellTypes_subset)



system.time(res_sub <- pbmcapply::pbmclapply(1:50, function(x) {
  
  l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), 
                      cellNames = colnames(exprsMat_subset), subset_test = F,  n = round(ncol(exprsMat_subset)/5)*4,  balance = F)
  table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)
}, mc.cores = 10))
gc()


res_seqwells <- do.call(cbind, res_sub)





################################################
########      10x Genomics       ###############
################################################

levinPBMCbenchmark_10x <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "10x Chromium (v3)"]

levinPBMCbenchmark_10x$batch <- droplevels(levinPBMCbenchmark_10x$Experiment)

logcounts(levinPBMCbenchmark_10x) <- as.matrix(logcounts(levinPBMCbenchmark_10x))
counts(levinPBMCbenchmark_10x) <- as.matrix(counts(levinPBMCbenchmark_10x))

levinPBMCbenchmark_10x <- levinPBMCbenchmark_10x[rowSums(counts(levinPBMCbenchmark_10x))!=0,]
levinPBMCbenchmark_10x 

exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_10x, "logcounts"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_10x$CellType)
table(cellTypes_subset)


set.seed(2019)



system.time(res_sub <- pbmcapply::pbmclapply(1:50, function(x) {
  
  l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), 
                      cellNames = colnames(exprsMat_subset), subset_test = F,  n = round(ncol(exprsMat_subset)/5)*4,  balance = F)
  table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)
}, mc.cores = 10))
gc()


res_10Xv3 <- do.call(cbind, res_sub)




################################################
########      10x Genomics (V2)      ###########
################################################

levinPBMCbenchmark_10xV2 <- levinPBMCbenchmark[, levinPBMCbenchmark$Method == "10x Chromium (v2)"]

levinPBMCbenchmark_10xV2$batch <- droplevels(levinPBMCbenchmark_10xV2$Experiment)

logcounts(levinPBMCbenchmark_10xV2) <- as.matrix(logcounts(levinPBMCbenchmark_10xV2))
counts(levinPBMCbenchmark_10xV2) <- as.matrix(counts(levinPBMCbenchmark_10xV2))

levinPBMCbenchmark_10xV2 <- levinPBMCbenchmark_10xV2[rowSums(counts(levinPBMCbenchmark_10xV2))!=0,]
levinPBMCbenchmark_10xV2 

exprsMat_subset <- as.matrix(assay(levinPBMCbenchmark_10xV2, "logcounts"))
exprsMat_subset <- as(exprsMat_subset, "sparseMatrix")
cellTypes_subset <- droplevels(levinPBMCbenchmark_10xV2$CellType)
table(cellTypes_subset)


system.time(res_sub <- pbmcapply::pbmclapply(1:50, function(x) {
  
  l <- subSamplingCal(exprsMat_subset, cellTypes_subset, geneNames = rownames(exprsMat_subset), 
                      cellNames = colnames(exprsMat_subset), subset_test = F,  n = round(ncol(exprsMat_subset)/5)*4,  balance = F)
  table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)
}, mc.cores = 10))
gc()


res_10Xv2 <- do.call(cbind, res_sub)






trainSelf_accuracy <- list(smartseq = res_smartseq,
                           celseq = res_celseq,
                           tenXv2 = res_10Xv2,
                           tenXv3 = res_10Xv3,
                           indrops = res_indrops,
                           dropseq = res_dropseq,
                           seqwells = res_seqwells
)


lapply(trainSelf_accuracy, rowMeans)

library(moon)
saveRDSwithDate("ensembleRes/trainSelf_accuracy", trainSelf_accuracy)

