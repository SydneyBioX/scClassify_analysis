### script to run DECENT

library(DECENT)
library(SingleCellExperiment)
library(scater)

##############################################################################
###################  Filter the doublet cells ################################
##############################################################################

sce_pbmc10k <- readRDS("/dskh/nobackup/yingxinl/ClassifyCell/data/sce_pbmc10k_v300.rds")
pbmc10k_doublet <- readRDS("results/pbmcResults/pbmc10k_doublet.rds")


sce_pbmc10k$doubletScores <- pbmc10k_doublet
sce_pbmc10k$doublet <- pbmc10k_doublet>0.5
sum(sce_pbmc10k$doublet)

table(sce_pbmc10k$doublet, sce_pbmc10k$Assigned_CellType)

plotTSNE(sce_pbmc10k, colour_by = "doubletScores")
plotTSNE(sce_pbmc10k, colour_by = "doublet")

sce_pbmc10k <- sce_pbmc10k[,!sce_pbmc10k$doublet]

# sce_pbmc10k <- runTSNE(sce_pbmc10k)
sce_pbmc10k$cellTypes <- sce_pbmc10k$Assigned_CellType
plotTSNE(sce_pbmc10k, colour_by = "cellTypes")

gc(reset = T)


##############################################################################
###################     Analysis               ###############################
##############################################################################





cellTypes_pbmc10k <- as.factor(sce_pbmc10k$cellTypes)

# levels(cellTypes_pbmc10k) <- c("B", "CD14+ Mono", "CD34+", 
#"Helper CD4 T", "NK", "CD8 T", 
# "Memory CD4 T", "Naive CD8 T", "Naive CD4 T", "Regulatory CD4 T")
table(cellTypes_pbmc10k)

sce_pbmc10k$cellTypes2 <- as.character(sce_pbmc10k$cellTypes)
sce_pbmc10k$cellTypes2[grep("CD4", sce_pbmc10k$cellTypes)] <- "CD4+ T"
sce_pbmc10k$cellTypes2[grep("CD8", sce_pbmc10k$cellTypes)] <- "CD8+ T"
plotTSNE(sce_pbmc10k, colour_by = "cellTypes2")


sce_pbmc10k$cellTypes3 <- as.character(sce_pbmc10k$cellTypes)
sce_pbmc10k$cellTypes3[grep("CD4+", sce_pbmc10k$cellTypes)] <- "T"
sce_pbmc10k$cellTypes3[grep("CD8+", sce_pbmc10k$cellTypes)] <- "T"
plotTSNE(sce_pbmc10k, colour_by = "cellTypes3")





##############################################################################
###################     DECENT                 ###############################
##############################################################################




# DECENT without spike-ins
de.table <- decent(data.obs = as.matrix(counts(sce_pbmc10k)),
                   X = ~as.factor(sce_pbmc10k$cellTypes), 
                   use.spikes = F,
                   CE.range = c(0.02, 0.1), # specify the range of the ranked random capture efficiency
                   dir = "DECENT/10x_pbmc10k_v3/",
                   parallel = TRUE,
                   n.cores = 8
)



decent.noDE <- readRDS("DECENT/10x_pbmc10k_v3/decent.noDE.rds")


getRho <- function(expr, obj) {
  # expr : expression (raw count) matrix. genes=row, cells=col
  # obj  : object from decent.noDE.rds
  # prob : thinning peobability (0-1)
  logit.rho <- obj$tau0[1] + obj$tau1[1]*log(expr + 0.1)
  rho       <- (1 + exp(-logit.rho))^-1
  rho
}


rhoMat <- getRho(as.matrix(counts(sce_pbmc10k)), decent.noDE)

getSubMat <- function(expr, rhoMat, prob, ncores = 20) {
  subMat <- pbmcapply::pbmclapply(1:ncol(expr), function(i) {
    VGAM::rbetabinom(nrow(expr),size = expr[,i],rho = rhoMat[, i], prob = prob)
  }, mc.cores = ncores)
  subMat <- do.call(cbind, subMat)
  return(subMat)
}

exprsMat_pbmc10k <- as.matrix(counts(sce_pbmc10k))
subMat <- getSubMat(exprsMat_pbmc10k, rhoMat, 0.1)








runDECENTsub_learningCurve <- function(exprsMat, cellTypes, rhoMat, prob, ncores = 5 
) {
  
  cat("Subset the matrix")
  subMat <- getSubMat(exprsMat, rhoMat, prob = prob, ncores = ncores)
  colnames(subMat) <- colnames(exprsMat)
  rownames(subMat) <- rownames(exprsMat)
  
  names(cellTypes) <- colnames(exprsMat)
  
  library(SingleCellExperiment)
  
  sce_sim <- SingleCellExperiment(assay = list(counts = subMat),
                                  colData = list(cellTypes = cellTypes))
  
  
  sce_sim <- scater::normalize(sce_sim)
  
  exprsMat_sim <- logcounts(sce_sim)
  
  cvIdx <- cvTools::cvFolds(ncol(sce_sim), K = 5)
  
  cellTypes_sim <- sce_sim$cellTypes
  cat("performing scClassify")
  exprsMat_sim <- as(exprsMat_sim, "sparseMatrix")
  
  n_list <- c(seq(20, 100, 20), seq(200, 2000, 200))
  
  res_sub <- list()
  for (i in 1:length(n_list)){ 
    print(paste("n=",n_list[i]))
    system.time(res_sub[[i]] <- pbmcapply::pbmclapply(1:5, function(x) {
      tryCatch({l <- subSamplingCal(exprsMat_sim, cellTypes_sim, geneNames = rownames(exprsMat_sim), 
                                    cellNames = colnames(exprsMat_sim),  
                                    subset_test = F,  n = n_list[[i]], balance = F, num_test = 1000)
      table(l$testRes$test$pearson_WKNN_limma$classifyRes)/length(l$testRes$test$pearson_WKNN_limma$classifyRes)}, 
      error = function(e){NULL})
    }, mc.cores = 5))
    print(do.call(cbind, res_sub[[i]]))
    print(paste("Median accuracy", mean(unlist(lapply(res_sub[[i]], "[[", "correct")))))
    gc()
  }
  
  
  names(res_sub) <- n_list
  res_sub <- lapply(res_sub, function(x) do.call(cbind, x))
  return(res_sub)
  

}







runDECENTsub <- function(exprsMat, cellTypes, rhoMat, prob, ncores = 5 
) {
  
  cat("Subset the matrix")
  subMat <- getSubMat(exprsMat, rhoMat, prob = prob, ncores = ncores)
  colnames(subMat) <- colnames(exprsMat)
  rownames(subMat) <- rownames(exprsMat)
  
  names(cellTypes) <- colnames(exprsMat)
  
  library(SingleCellExperiment)
  
  sce_sim <- SingleCellExperiment(assay = list(counts = subMat),
                                  colData = list(cellTypes = cellTypes))
  
  
  sce_sim <- scater::normalize(sce_sim)
  
  exprsMat_sim <- logcounts(sce_sim)
  
  cvIdx <- cvTools::cvFolds(ncol(sce_sim), K = 5)
  
  cat("performing scClassify")
  
  
  classRes_list <- pbmcapply::pbmclapply(1:5, function(cv){
    trainIdx <- cvIdx$subsets[cvIdx$which != cv]
    # trainIdx <- sample(ncol(exprsMat_sim), round(ncol(exprsMat_sim)) * 0.8)
    table(sce_sim$cellTypes[trainIdx])
    
    
    
    testRes <- scClassify(exprsMat_train = exprsMat_sim[, trainIdx],
                          cellTypes_train = sce_sim$cellTypes[trainIdx],
                          exprsMat_test = list(test = exprsMat_sim[, -trainIdx]),
                          cellTypes_test = list(test = sce_sim$cellTypes[-trainIdx]),
                          tree = "HOPACH",
                          algorithm = "WKNN",
                          k = 5,
                          selectFeatures = c("limma"),
                          similarity = c("pearson"),
                          verbose = F)
    
    
    classRes <- testRes$testRes$test$pearson_WKNN_limma$classifyRes
    table(classRes)/length(classRes)
  }, mc.cores = 5)
  
  res <- do.call(cbind, classRes_list)
  print(res)
  return(res)
}





#### Level 2

rho_res_level2 <- list()
for (i in 1:length(rho_list)) {
  print(paste("===========", rho_list[i], "========"))
  rho_res_level2[[i]] <- pbmcapply::pbmclapply(1:20, function(x) 
    runDECENTsub(exprsMat_pbmc10k, sce_pbmc10k$cellTypes2, rhoMat, prob = rho_list[i]), mc.cores = 4)
}


rho_res_level2 <- lapply(rho_res_level2, function(x) do.call(cbind, x))

names(rho_res_level2) <- rho_list


saveRDS(rho_res_level2, "results/decent_rho_res_level2_pbmc10k.rds")


library(reshape2)
get_rho_res <- melt(rho_res_level2)

ggplot(get_rho_res[get_rho_res$Var1 == "correct", ], aes(x = L1, y = value)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Rho") +
  ylab("Accuracy Rate")



#### Level 1

rho_res_level1 <- list()
for (i in 1:length(rho_list)) {
  print(paste("===========", rho_list[i], "========"))
  rho_res_level1[[i]] <- pbmcapply::pbmclapply(1:20, function(x) 
    runDECENTsub(exprsMat_pbmc10k, sce_pbmc10k$cellTypes3, rhoMat, prob = rho_list[i]), mc.cores = 4)
}


rho_res_level1 <- lapply(rho_res_level1, function(x) do.call(cbind, x))

names(rho_res_level1) <- rho_list

saveRDS(rho_res_level1, "results/decent_rho_res_level1_pbmc10k.rds")

library(reshape2)
get_rho_res <- melt(rho_res_level1)

ggplot(get_rho_res[get_rho_res$Var1 == "correct", ], aes(x = L1, y = value)) +
  geom_boxplot() +
  theme_bw() +
  xlab("Rho") +
  ylab("Accuracy Rate")




get_rho_res_level1 <- melt(rho_res_level1)
get_rho_res_level1$level <- "Level1"
get_rho_res_level2 <- melt(rho_res_level2)
get_rho_res_level2$level <- "Level2"



rho_list <- c(seq(0.1, 1, 0.1))

start <- Sys.time()
res <- runDECENTsub_learningCurve(exprsMat_pbmc10k[, 1:2200], sce_pbmc10k$cellTypes3[1:2200], rhoMat[, 1:2200], prob = 0.5)
end <- Sys.time()
start - end


rho_list <- c(0.1, 0.2, 0.5, 0.8, 1)


rho_res_level2 <- list()
for (i in 1:length(rho_list)) {
  print(paste("===========", rho_list[i], "========"))
  rho_res_level2[[i]] <- pbmcapply::pbmclapply(1:20, function(x) 
    runDECENTsub_learningCurve(exprsMat_pbmc10k, sce_pbmc10k$cellTypes2, rhoMat, prob = rho_list[i]), 
    mc.cores = 5)
  print(rho_res_level2[[i]][[1]])
}




saveRDS(rho_res_level2,
  file = "results/res_decent_pbmc10k_level2_learningCurveRes.rds")




rho_res_level1 <- list()
for (i in 1:length(rho_list)) {
  print(paste("===========", rho_list[i], "========"))
  rho_res_level1[[i]] <- pbmcapply::pbmclapply(1:20, function(x) 
    runDECENTsub_learningCurve(exprsMat_pbmc10k, sce_pbmc10k$cellTypes3, rhoMat, prob = rho_list[i]), 
    mc.cores = 5)
  print(rho_res_level1[[i]][[1]])
}


saveRDS(rho_res_level1,
  file = "results/res_decent_pbmc10k_level1_learningCurveRes.rds")

