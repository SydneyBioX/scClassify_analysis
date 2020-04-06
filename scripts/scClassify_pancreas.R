
set.seed(1)


# Package


library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(CHETAH)
library(pheatmap)

library(limma)
library(hopach)
library(ggraph)
library(igraph)
library(mixtools)
library(pbmcapply)
library(proxy)
library(minpack.lm)

rscript <- list.files("scClassify/")
sapply(rscript, function(x){
  source(paste("scClassify/", x, sep = ""))
})



saveRDSwithDate <- function(fileName, ...) {
  today <- gsub("-","", Sys.Date())
  fileName <- paste(paste(fileName, today, sep = "_"), ".rds", sep = "")
  saveRDS(..., file = fileName)
}


# Data 

load("data/pancreasSeven.RData")


baron <- baron[,grep("human", colnames(baron))]

baron
muraro1
muraro2
segerstolpe
lawlor
xin
wang

xin$cellTypes[xin$cellTypes=="PP"] <- "gamma"

table(baron$cellTypes)
table(muraro1$cellTypes)
table(muraro2$cellTypes)
table(segerstolpe$cellTypes)
table(lawlor$cellTypes)
table(xin$cellTypes)
table(wang$cellTypes)

muraro2$cellTypes[muraro2$cellTypes == "mesenchymal"] <- "stellate"

wang$cellTypes[wang$cellTypes == "mesenchymal"] <- "stellate"

segerstolpe <- segerstolpe[,!segerstolpe$cellTypes%in%c("co-expression","unclassified endocrine")]






exprsMat_wang <- logcounts(wang)
exprsMat_xin <- logcounts(xin)
exprsMat_lawlor <- logcounts(lawlor)
exprsMat_segerstolpe <- logcounts(segerstolpe)
exprsMat_muraro <- logcounts(muraro2)
exprsMat_baron <- logcounts(baron)




system.time(wang_scClassify_res <- scClassify(exprsMat_train = exprsMat_wang,
                                                       cellTypes_train = wang$cellTypes,
                                                       exprsMat_test = list(xin = exprsMat_xin,
                                                                            lawlor = exprsMat_lawlor,
                                                                            segerstolpe = exprsMat_segerstolpe,
                                                                            baron = exprsMat_baron,
                                                                            muraro = exprsMat_muraro),
                                                       cellTypes_test = list(xin = xin$cellTypes,
                                                                             lawlor = lawlor$cellTypes,
                                                                             segerstolpe = segerstolpe$cellTypes,
                                                                             baron = baron$cellTypes,
                                                                             muraro = muraro2$cellTypes),
                                                       similarity = c("pearson"),
                                                       selectFeatures = c("limma"),
                                                       algorithm = c("WKNN"),
                                                       parallel = F,
                                                       ncores = 10))

do.call(rbind, lapply(wang_scClassify_res$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

system.time(lawlor_scClassify_res <- scClassify(exprsMat_train = exprsMat_lawlor,
                                                         cellTypes_train = lawlor$cellTypes,
                                                         exprsMat_test = list(xin = exprsMat_xin,
                                                                              wang = exprsMat_wang,
                                                                              segerstolpe = exprsMat_segerstolpe,
                                                                              baron = exprsMat_baron,
                                                                              muraro = exprsMat_muraro),
                                                         cellTypes_test = list(xin = xin$cellTypes,
                                                                               wang = wang$cellTypes,
                                                                               segerstolpe = segerstolpe$cellTypes,
                                                                               baron = baron$cellTypes,
                                                                               muraro = muraro2$cellTypes),
                                                similarity = c("pearson"),
                                                selectFeatures = c("limma"),
                                                algorithm = c("WKNN"),
                                                parallel = F,
                                                         ncores = 10))
do.call(rbind, lapply(lawlor_scClassify_res$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

segerstolpe_scClassify_res <- scClassify(exprsMat_train = exprsMat_segerstolpe,
                                                  cellTypes_train = segerstolpe$cellTypes,
                                                  exprsMat_test = list(xin = exprsMat_xin,
                                                                       wang = exprsMat_wang,
                                                                       lawlor = exprsMat_lawlor,
                                                                       baron = exprsMat_baron,
                                                                       muraro = exprsMat_muraro),
                                                  cellTypes_test = list(xin = xin$cellTypes,
                                                                        wang = wang$cellTypes,
                                                                        lawlor = lawlor$cellTypes,
                                                                        baron = baron$cellTypes,
                                                                        muraro = muraro2$cellTypes),
                                                  similarity = c("pearson"),
                                                  selectFeatures = c("limma"),
                                                  algorithm = c("WKNN"),
                                                  parallel = F,
                                                  ncores = 10)

do.call(rbind, lapply(segerstolpe_scClassify_res$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

xin_scClassify_res <- scClassify(exprsMat_train = exprsMat_xin,
                                          cellTypes_train = xin$cellTypes,
                                          exprsMat_test = list(segerstolpe = exprsMat_segerstolpe,
                                                               wang = exprsMat_wang,
                                                               lawlor = exprsMat_lawlor,
                                                               baron = exprsMat_baron,
                                                               muraro = exprsMat_muraro),
                                          cellTypes_test = list(segerstolpe = segerstolpe$cellTypes,
                                                                wang = wang$cellTypes,
                                                                lawlor = lawlor$cellTypes,
                                                                baron = baron$cellTypes,
                                                                muraro = muraro2$cellTypes),
                                          similarity = c("pearson"),
                                          selectFeatures = c("limma"),
                                          algorithm = c("WKNN"),
                                          parallel = F,
                                          ncores = 10)

do.call(rbind, lapply(xin_scClassify_res$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))


muraro_scClassify_res <- scClassify(exprsMat_train = exprsMat_muraro,
                                             cellTypes_train = muraro2$cellTypes,
                                             exprsMat_test = list(segerstolpe = exprsMat_segerstolpe,
                                                                  wang = exprsMat_wang,
                                                                  lawlor = exprsMat_lawlor,
                                                                  baron = exprsMat_baron,
                                                                  xin = exprsMat_xin),
                                             cellTypes_test = list(segerstolpe = segerstolpe$cellTypes,
                                                                   wang = wang$cellTypes,
                                                                   lawlor = lawlor$cellTypes,
                                                                   baron = baron$cellTypes,
                                                                   xin = xin$cellTypes),
                                             similarity = c("pearson"),
                                             selectFeatures = c("limma"),
                                             algorithm = c("WKNN"),
                                             parallel = F,
                                             ncores = 10)
do.call(rbind, lapply(muraro_scClassify_res$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

baron_scClassify_res <- scClassify(exprsMat_train = exprsMat_baron,
                                            cellTypes_train = baron$cellTypes,
                                            exprsMat_test = list(segerstolpe = exprsMat_segerstolpe,
                                                                 wang = exprsMat_wang,
                                                                 lawlor = exprsMat_lawlor,
                                                                 xin = exprsMat_xin,
                                                                 muraro = exprsMat_muraro),
                                            cellTypes_test = list(segerstolpe = segerstolpe$cellTypes,
                                                                  wang = wang$cellTypes,
                                                                  lawlor = lawlor$cellTypes,
                                                                  xin = xin$cellTypes,
                                                                  muraro = muraro2$cellTypes),
                                            similarity = c("pearson"),
                                            selectFeatures = c("limma"),
                                            algorithm = c("WKNN"),
                                            parallel = F,
                                            ncores = 10)


do.call(rbind, lapply(baron_scClassify_res$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))


scClassifyRes <- rbind(
  do.call(rbind, lapply(baron_scClassify_res$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(muraro_scClassify_res$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(segerstolpe_scClassify_res$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(lawlor_scClassify_res$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(wang_scClassify_res$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(xin_scClassify_res$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))
  
)

rownames(scClassifyRes) <- paste(c(rep("baron", 5),
                                   rep("muraro2", 5),
                                   rep("segerstolpe", 5), 
                                   rep("lawlor", 5),
                                   rep("wang", 5),
                                   rep("xin", 5)), rownames(scClassifyRes), sep = "_")




system.time(wang_scClassify_res_prob70 <- scClassify(exprsMat_train = exprsMat_wang,
                                              cellTypes_train = wang$cellTypes,
                                              exprsMat_test = list(xin = exprsMat_xin,
                                                                   lawlor = exprsMat_lawlor,
                                                                   segerstolpe = exprsMat_segerstolpe,
                                                                   baron = exprsMat_baron,
                                                                   muraro = exprsMat_muraro),
                                              cellTypes_test = list(xin = xin$cellTypes,
                                                                    lawlor = lawlor$cellTypes,
                                                                    segerstolpe = segerstolpe$cellTypes,
                                                                    baron = baron$cellTypes,
                                                                    muraro = muraro2$cellTypes),
                                              similarity = c("pearson"),
                                              selectFeatures = c("limma"),
                                              algorithm = c("WKNN"),
                                              parallel = F,
                                              prob_threshold = 0.7,
                                              ncores = 10))

do.call(rbind, lapply(wang_scClassify_res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

system.time(lawlor_scClassify_res_prob70 <- scClassify(exprsMat_train = exprsMat_lawlor,
                                                cellTypes_train = lawlor$cellTypes,
                                                exprsMat_test = list(xin = exprsMat_xin,
                                                                     wang = exprsMat_wang,
                                                                     segerstolpe = exprsMat_segerstolpe,
                                                                     baron = exprsMat_baron,
                                                                     muraro = exprsMat_muraro),
                                                cellTypes_test = list(xin = xin$cellTypes,
                                                                      wang = wang$cellTypes,
                                                                      segerstolpe = segerstolpe$cellTypes,
                                                                      baron = baron$cellTypes,
                                                                      muraro = muraro2$cellTypes),
                                                similarity = c("pearson"),
                                                selectFeatures = c("limma"),
                                                algorithm = c("WKNN"),
                                                parallel = F,
                                                prob_threshold = 0.7,
                                                ncores = 10))
do.call(rbind, lapply(lawlor_scClassify_res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

segerstolpe_scClassify_res_prob70 <- scClassify(exprsMat_train = exprsMat_segerstolpe,
                                         cellTypes_train = segerstolpe$cellTypes,
                                         exprsMat_test = list(xin = exprsMat_xin,
                                                              wang = exprsMat_wang,
                                                              lawlor = exprsMat_lawlor,
                                                              baron = exprsMat_baron,
                                                              muraro = exprsMat_muraro),
                                         cellTypes_test = list(xin = xin$cellTypes,
                                                               wang = wang$cellTypes,
                                                               lawlor = lawlor$cellTypes,
                                                               baron = baron$cellTypes,
                                                               muraro = muraro2$cellTypes),
                                         similarity = c("pearson"),
                                         selectFeatures = c("limma"),
                                         algorithm = c("WKNN"),
                                         parallel = F,
                                         prob_threshold = 0.7,
                                         ncores = 10)

do.call(rbind, lapply(segerstolpe_scClassify_res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

xin_scClassify_res_prob70 <- scClassify(exprsMat_train = exprsMat_xin,
                                 cellTypes_train = xin$cellTypes,
                                 exprsMat_test = list(segerstolpe = exprsMat_segerstolpe,
                                                      wang = exprsMat_wang,
                                                      lawlor = exprsMat_lawlor,
                                                      baron = exprsMat_baron,
                                                      muraro = exprsMat_muraro),
                                 cellTypes_test = list(segerstolpe = segerstolpe$cellTypes,
                                                       wang = wang$cellTypes,
                                                       lawlor = lawlor$cellTypes,
                                                       baron = baron$cellTypes,
                                                       muraro = muraro2$cellTypes),
                                 similarity = c("pearson"),
                                 selectFeatures = c("limma"),
                                 algorithm = c("WKNN"),
                                 parallel = F,
                                 prob_threshold = 0.7,
                                 ncores = 10)

do.call(rbind, lapply(xin_scClassify_res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))


muraro_scClassify_res_prob70 <- scClassify(exprsMat_train = exprsMat_muraro,
                                    cellTypes_train = muraro2$cellTypes,
                                    exprsMat_test = list(segerstolpe = exprsMat_segerstolpe,
                                                         wang = exprsMat_wang,
                                                         lawlor = exprsMat_lawlor,
                                                         baron = exprsMat_baron,
                                                         xin = exprsMat_xin),
                                    cellTypes_test = list(segerstolpe = segerstolpe$cellTypes,
                                                          wang = wang$cellTypes,
                                                          lawlor = lawlor$cellTypes,
                                                          baron = baron$cellTypes,
                                                          xin = xin$cellTypes),
                                    similarity = c("pearson"),
                                    selectFeatures = c("limma"),
                                    algorithm = c("WKNN"),
                                    parallel = F,
                                    prob_threshold = 0.7,
                                    ncores = 10)
do.call(rbind, lapply(muraro_scClassify_res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))

baron_scClassify_res_prob70 <- scClassify(exprsMat_train = exprsMat_baron,
                                   cellTypes_train = baron$cellTypes,
                                   exprsMat_test = list(segerstolpe = exprsMat_segerstolpe,
                                                        wang = exprsMat_wang,
                                                        lawlor = exprsMat_lawlor,
                                                        xin = exprsMat_xin,
                                                        muraro = exprsMat_muraro),
                                   cellTypes_test = list(segerstolpe = segerstolpe$cellTypes,
                                                         wang = wang$cellTypes,
                                                         lawlor = lawlor$cellTypes,
                                                         xin = xin$cellTypes,
                                                         muraro = muraro2$cellTypes),
                                   similarity = c("pearson"),
                                   selectFeatures = c("limma"),
                                   algorithm = c("WKNN"),
                                   parallel = F,
                                   prob_threshold = 0.7,
                                   ncores = 10)


do.call(rbind, lapply(baron_scClassify_res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))


scClassifyRes_prob70 <- rbind(
  do.call(rbind, lapply(baron_scClassify_res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(muraro_scClassify_res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(segerstolpe_scClassify_res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(lawlor_scClassify_res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(wang_scClassify_res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes))),
  do.call(rbind, lapply(xin_scClassify_res_prob70$testRes, function(x) table(x$pearson_WKNN_limma$classifyRes)/length(x$pearson_WKNN_limma$classifyRes)))
  
)

rownames(scClassifyRes_prob70) <- paste(c(rep("baron", 5),
                                   rep("muraro2", 5),
                                   rep("segerstolpe", 5), 
                                   rep("lawlor", 5),
                                   rep("wang", 5),
                                   rep("xin", 5)), rownames(scClassifyRes_prob70), sep = "_")


saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/results/pancreasResults/pancreasSix_scClassifyRes_prob70", scClassifyRes_prob70)





###################################################################################################
####################################### ENSEMBLE ##################################################
###################################################################################################




system.time(wang_scClassify_ensemble_res <- scClassify(exprsMat_train = exprsMat_wang,
                                                       cellTypes_train = wang$cellTypes,
                                                       exprsMat_test = list(xin = exprsMat_xin,
                                                                            lawlor = exprsMat_lawlor,
                                                                            segerstolpe = exprsMat_segerstolpe,
                                                                            baron = exprsMat_baron,
                                                                            muraro = exprsMat_muraro),
                                                       cellTypes_test = list(xin = xin$cellTypes,
                                                                             lawlor = lawlor$cellTypes,
                                                                             segerstolpe = segerstolpe$cellTypes,
                                                                             baron = baron$cellTypes,
                                                                             muraro = muraro2$cellTypes),
                                                       similarity = c("pearson",  "spearman",
                                                                      "cosine", "jaccard", "kendall",
                                                                      "weighted_rank"),
                                                       selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                       algorithm = c("WKNN", "KNN"),
                                                       parallel = T,
                                                       prob_threshold = 0.7,
                                                       ncores = 20))

saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/wang_scClassify_ensemble_res", wang_scClassify_ensemble_res)

system.time(lawlor_scClassify_ensemble_res <- scClassify(exprsMat_train = exprsMat_lawlor,
                                                         cellTypes_train = lawlor$cellTypes,
                                                         exprsMat_test = list(xin = exprsMat_xin,
                                                                              wang = exprsMat_wang,
                                                                              segerstolpe = exprsMat_segerstolpe,
                                                                              baron = exprsMat_baron,
                                                                              muraro = exprsMat_muraro),
                                                         cellTypes_test = list(xin = xin$cellTypes,
                                                                               wang = wang$cellTypes,
                                                                               segerstolpe = segerstolpe$cellTypes,
                                                                               baron = baron$cellTypes,
                                                                               muraro = muraro2$cellTypes),
                                                         similarity = c("pearson",  "spearman",
                                                                        "cosine", "jaccard", "kendall",
                                                                        "weighted_rank"),
                                                         selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                         algorithm = c("WKNN", "KNN"),
                                                         parallel = T,
                                                         prob_threshold = 0.7,
                                                         ncores = 20))

saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/lawlor_scClassify_ensemble_res", lawlor_scClassify_ensemble_res)

segerstolpe_scClassify_ensemble_res <- scClassify(exprsMat_train = exprsMat_segerstolpe,
                                                  cellTypes_train = segerstolpe$cellTypes,
                                                  exprsMat_test = list(xin = exprsMat_xin,
                                                                       wang = exprsMat_wang,
                                                                       lawlor = exprsMat_lawlor,
                                                                       baron = exprsMat_baron,
                                                                       muraro = exprsMat_muraro),
                                                  cellTypes_test = list(xin = xin$cellTypes,
                                                                        wang = wang$cellTypes,
                                                                        lawlor = lawlor$cellTypes,
                                                                        baron = baron$cellTypes,
                                                                        muraro = muraro2$cellTypes),
                                                  similarity = c("pearson",  "spearman",
                                                                 "cosine", "jaccard", "kendall",
                                                                 "weighted_rank"),
                                                  selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                  algorithm = c("WKNN", "KNN"),
                                                  parallel = T,
                                                  prob_threshold = 0.7,
                                                  ncores = 20)


saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/segerstolpe_scClassify_ensemble_res", segerstolpe_scClassify_ensemble_res)

xin_scClassify_ensemble_res <- scClassify(exprsMat_train = exprsMat_xin,
                                          cellTypes_train = xin$cellTypes,
                                          exprsMat_test = list(segerstolpe = exprsMat_segerstolpe,
                                                               wang = exprsMat_wang,
                                                               lawlor = exprsMat_lawlor,
                                                               baron = exprsMat_baron,
                                                               muraro = exprsMat_muraro),
                                          cellTypes_test = list(segerstolpe = segerstolpe$cellTypes,
                                                                wang = wang$cellTypes,
                                                                lawlor = lawlor$cellTypes,
                                                                baron = baron$cellTypes,
                                                                muraro = muraro2$cellTypes),
                                          similarity = c("pearson",  "spearman",
                                                         "cosine", "jaccard", "kendall",
                                                         "weighted_rank"),
                                          selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                          algorithm = c("WKNN", "KNN"),
                                          parallel = T,
                                          prob_threshold = 0.7,
                                          ncores = 20)


saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/xin_scClassify_ensemble_res", xin_scClassify_ensemble_res)

muraro_scClassify_ensemble_res <- scClassify(exprsMat_train = exprsMat_muraro,
                                             cellTypes_train = muraro2$cellTypes,
                                             exprsMat_test = list(segerstolpe = exprsMat_segerstolpe,
                                                                  wang = exprsMat_wang,
                                                                  lawlor = exprsMat_lawlor,
                                                                  baron = exprsMat_baron,
                                                                  xin = exprsMat_xin),
                                             cellTypes_test = list(segerstolpe = segerstolpe$cellTypes,
                                                                   wang = wang$cellTypes,
                                                                   lawlor = lawlor$cellTypes,
                                                                   baron = baron$cellTypes,
                                                                   xin = xin$cellTypes),
                                             similarity = c("pearson",  "spearman",
                                                            "cosine", "jaccard", "kendall",
                                                            "weighted_rank"),
                                             selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                             algorithm = c("WKNN", "KNN"),
                                             parallel = T,
                                             prob_threshold = 0.7,
                                             ncores = 20)

saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/muraro_scClassify_ensemble_res", muraro_scClassify_ensemble_res)

baron_scClassify_ensemble_res <- scClassify(exprsMat_train = exprsMat_baron,
                                            cellTypes_train = baron$cellTypes,
                                            exprsMat_test = list(segerstolpe = exprsMat_segerstolpe,
                                                                 wang = exprsMat_wang,
                                                                 lawlor = exprsMat_lawlor,
                                                                 xin = exprsMat_xin,
                                                                 muraro = exprsMat_muraro),
                                            cellTypes_test = list(segerstolpe = segerstolpe$cellTypes,
                                                                  wang = wang$cellTypes,
                                                                  lawlor = lawlor$cellTypes,
                                                                  xin = xin$cellTypes,
                                                                  muraro = muraro2$cellTypes),
                                            similarity = c("pearson",  "spearman",
                                                           "cosine", "jaccard", "kendall",
                                                           "weighted_rank"),
                                            selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                            algorithm = c("WKNN", "KNN"),
                                            parallel = T,
                                            prob_threshold = 0.7,
                                            ncores = 20)


saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/baron_scClassify_ensemble_res", baron_scClassify_ensemble_res)

  
  

library(cvTools)
scClassify_trainAccuracy <- function(exprsMat_train,
                                     cellTypes_train,
                                     cvFolds = 5,
                                     parallel = T,
                                     ncores = 10
) {
  
  cv <- cvTools::cvFolds(ncol(exprsMat_train), K = cvFolds)
  
  cvRes <- list()
  for (cvIdx in 1:cvFolds) {
    trainIdx <- cv$subsets[cv$which != cvIdx]
    testIdx <-  cv$subsets[cv$which == cvIdx]
    cvRes[[cvIdx]] <- scClassify(exprsMat_train = exprsMat_train[,trainIdx],
                                 cellTypes_train = cellTypes_train[trainIdx],
                                 exprsMat_test = exprsMat_train[,testIdx],
                                 cellTypes_test = cellTypes_train[testIdx],
                                 similarity = c("pearson",  "spearman",
                                                "cosine", "jaccard", "kendall",
                                                "weighted_rank","manhattan"),
                                 selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                 algorithm = c("WKNN", "KNN"),
                                 parallel = parallel,
                                 ncores = ncores)
  }
  return(cvRes)
  
}

getEnsembleError <- function(res){
  ensembleErr <- do.call(rbind, lapply(res, function(x) table(x$classifyRes)/length(x$classifyRes)))
  return(ensembleErr)
}





system.time(wang_scClassify_ensemble_res_itself <- scClassify(exprsMat_train = exprsMat_wang,
                                                       cellTypes_train = wang$cellTypes,
                                                       exprsMat_test = list(wang = exprsMat_wang
                                                                            ),
                                                       cellTypes_test = list(wang = wang$cellTypes),
                                                       similarity = c("pearson",  "spearman",
                                                                      "cosine", "jaccard", "kendall",
                                                                      "weighted_rank"),
                                                       selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                       algorithm = c("WKNN", "KNN"),
                                                       parallel = T,
                                                       prob_threshold = 0.7,
                                                       ncores = 10))

saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/wang_scClassify_ensemble_res_itself", wang_scClassify_ensemble_res_itself)

system.time(lawlor_scClassify_ensemble_res_itself <- scClassify(exprsMat_train = exprsMat_lawlor,
                                                         cellTypes_train = lawlor$cellTypes,
                                                         exprsMat_test = list(lawlor = exprsMat_lawlor),
                                                         cellTypes_test = list(lawlor = lawlor$cellTypes),
                                                         similarity = c("pearson",  "spearman",
                                                                        "cosine", "jaccard", "kendall",
                                                                        "weighted_rank"),
                                                         selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                         algorithm = c("WKNN", "KNN"),
                                                         parallel = T,
                                                         prob_threshold = 0.7,
                                                         ncores = 10))

saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/lawlor_scClassify_ensemble_res_itself", lawlor_scClassify_ensemble_res_itself)

segerstolpe_scClassify_ensemble_res_itself <- scClassify(exprsMat_train = exprsMat_segerstolpe,
                                                  cellTypes_train = segerstolpe$cellTypes,
                                                  exprsMat_test = list(segerstolpe = exprsMat_segerstolpe),
                                                  cellTypes_test = list(segerstolpe = segerstolpe$cellTypes),
                                                  similarity = c("pearson",  "spearman",
                                                                 "cosine", "jaccard", "kendall",
                                                                 "weighted_rank"),
                                                  selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                                  algorithm = c("WKNN", "KNN"),
                                                  parallel = T,
                                                  prob_threshold = 0.7,
                                                  ncores = 10)


saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/segerstolpe_scClassify_ensemble_res_itself", segerstolpe_scClassify_ensemble_res_itself)

xin_scClassify_ensemble_res_itself <- scClassify(exprsMat_train = exprsMat_xin,
                                          cellTypes_train = xin$cellTypes,
                                          exprsMat_test = list(xin = exprsMat_xin),
                                          cellTypes_test = list(xin = xin$cellTypes),
                                          similarity = c("pearson",  "spearman",
                                                         "cosine", "jaccard", "kendall",
                                                         "weighted_rank"),
                                          selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                          algorithm = c("WKNN", "KNN"),
                                          parallel = T,
                                          prob_threshold = 0.7,
                                          ncores = 15)


saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/xin_scClassify_ensemble_res_itself", xin_scClassify_ensemble_res_itself)

muraro_scClassify_ensemble_res_itself <- scClassify(exprsMat_train = exprsMat_muraro,
                                             cellTypes_train = muraro2$cellTypes,
                                             exprsMat_test = list(muraro = exprsMat_muraro),
                                             cellTypes_test = list(muraro = muraro2$cellTypes),
                                             similarity = c("pearson",  "spearman",
                                                            "cosine", "jaccard", "kendall",
                                                            "weighted_rank"),
                                             selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                             algorithm = c("WKNN", "KNN"),
                                             parallel = T,
                                             prob_threshold = 0.7,
                                             ncores = 15)

saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/muraro_scClassify_ensemble_res_itself", muraro_scClassify_ensemble_res_itself)

baron_scClassify_ensemble_res_itself <- scClassify(exprsMat_train = exprsMat_baron,
                                            cellTypes_train = baron$cellTypes,
                                            exprsMat_test = list(baron = exprsMat_baron),
                                            cellTypes_test = list(baron = baron$cellTypes),
                                            similarity = c("pearson",  "spearman",
                                                           "cosine", "jaccard", "kendall",
                                                           "weighted_rank"),
                                            selectFeatures = c("limma", "DV", "DD", "chisq", "BI"),
                                            algorithm = c("WKNN", "KNN"),
                                            parallel = T,
                                            prob_threshold = 0.7,
                                            ncores = 10)


saveRDSwithDate("/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/baron_scClassify_ensemble_res_itself", baron_scClassify_ensemble_res_itself)







lawlorWeights <- lapply(1:5, scClassify_trainAccuracy(exprsMat_train = exprsMat_lawlor,
                                                       cellTypes_train = lawlor$cellTypes))

saveRDS(lawlorWeights, file = "/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/lawlor_scClassify_ensemble_weights20190615.Rds")

wangWeights <- lapply(1:20, scClassify_trainAccuracy(exprsMat_train = exprsMat_wang,
                                                     cellTypes_train = wang$cellTypes))

saveRDS(wangWeights, file = "/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/wang_scClassify_ensemble_weights20190615.Rds")

xinWeights <- lapply(1:20, scClassify_trainAccuracy(exprsMat_train = exprsMat_xin,
                                                    cellTypes_train = xin$cellTypes))
saveRDS(xinWeights, file = "/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/xin_scClassify_ensemble_weights20190615.Rds")


segerstolpeWeights <-  lapply(1:20, scClassify_trainAccuracy(exprsMat_train = exprsMat_segerstolpe,
                                                             cellTypes_train = segerstolpe$cellTypes))
saveRDS(segerstolpeWeights, file = "/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/segerstolpe_scClassify_ensemble_weights20190615.Rds")



muraroWeights <-  lapply(1:20, scClassify_trainAccuracy(exprsMat_train = exprsMat_muraro,
                                                        cellTypes_train = muraro2$cellTypes))
saveRDS(muraroWeights, file = "/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/muraro_scClassify_ensemble_weights20190615.Rds")


baronWeights <-  lapply(1:20, scClassify_trainAccuracy(exprsMat_train = exprsMat_baron,
                                                       cellTypes_train = baron$cellTypes))
saveRDS(baronWeights, file = "/verona/nobackup/yingxinl/ClassifyCell/ensembleRes/baron_scClassify_ensemble_weights20190615.Rds")
