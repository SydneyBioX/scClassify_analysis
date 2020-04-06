rm(list = ls())


library(SingleCellExperiment)
library(scater)
library(scran)
library(ggplot2)
library(gridExtra)
library(ggthemes)
library(CHETAH)
library(pheatmap)
library(scdney)
source("/albona/nobackup/yingxinl/ClassifyCell/scClassify/featureSelection.R")

rscript <- list.files("scClassify/")
sapply(rscript, function(x){
  source(paste("scClassify/", x, sep = ""))
})


library(colortools)
cellTypes_col <- toupper(c(tableau_color_pal("Tableau 20")(20)[-c(13:14)],
                           tableau_color_pal("Classic 20")(20)[-c(15:16)]))

some_colors <- setColors(tableau_color_pal("Classic 20")(1), 70)

pizza(some_colors)


tabulaMuris_facs <- readRDS("data/tabulaMuris_facs.rds")

tabulaMuris_dropseq <- readRDS("data/tabulaMuris_dropseq.rds")



tabulaMuris_facs
tabulaMuris_dropseq

tabulaMuris_facs <- tabulaMuris_facs[rowSums(counts(tabulaMuris_facs))!= 0,]
tabulaMuris_dropseq <- tabulaMuris_dropseq[rowSums(counts(tabulaMuris_dropseq)) != 0,]


tabulaMuris_facs
tabulaMuris_dropseq


system.time(hvg <- featureSelection(as.matrix(logcounts(tabulaMuris_facs)), tabulaMuris_facs$cell_ontology_class,
feature = "limma", topN = 100))
length(na.omit(hvg))

trainFACS <- scClassify(exprsMat_train = as.matrix(logcounts(tabulaMuris_facs)[na.omit(hvg), ]),
                        cellTypes_train = tabulaMuris_facs$cell_ontology_class,
                        exprsMat_test = list(
                          tm = as.matrix(logcounts(tabulaMuris_dropseq)),
                          tm_facs = as.matrix(logcounts(tabulaMuris_facs))),
                        cellTypes_test = list(
                          tm = tabulaMuris_dropseq$cell_ontology_class,
                          tm_facs = tabulaMuris_facs$cell_ontology_class
                        ),
                        similarity = c("pearson"),
                        selectFeatures = c("limma"),
                        algorithm = c("WKNN"),
                        prob_threshold = 0.7,
                        hopach_kmax = 9,
                        verbose = T,
                        k = 10,
                        topN = 50)



saveRDS(trainFACS, "results/facs_trainRes")


### lung

rm(list = ls())

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
library(scdney)


tabulaMuris_dropseq_lung <- readRDS("data/tabulaMuris_dropseq_lung.rds")
tabulaMuris_facs_lung <- readRDS("data/tabulaMuris_facs_lung.rds")
cohen_lung <- readRDS("data/GSE119228_lung_sce.rds")



cohen_lung$cellTypes <- as.character(cohen_lung$cellTypes)
cohen_lung$cellTypes[cohen_lung$cellTypes == "Mast"] <- "Mast cell"
cohen_lung$cellTypes[cohen_lung$cellTypes == "Fibro"] <- "stromal cell"
cohen_lung$cellTypes[cohen_lung$cellTypes == "Ciliated"] <- "ciliated cell"
cohen_lung$cellTypes[cohen_lung$cellTypes == "NK"] <- "natural killer cell"
cohen_lung$cellTypes[cohen_lung$cellTypes == "Clara"] <- "Clara cell"
cohen_lung$cellTypes[cohen_lung$cellTypes == "Neut"] <- "Neutrophils"
cohen_lung$cellTypes[cohen_lung$cellTypes == "Smooth"] <- "Smooth muscle cells"
cohen_lung$cellTypes[cohen_lung$cellTypes == "T"] <- "T cell"
cohen_lung$cellTypes[cohen_lung$cellTypes == "B"] <- "B cell"
cohen_lung$cellTypes[cohen_lung$cellTypes == "DC"] <- "dendritic cell"
cohen_lung$cellTypes[cohen_lung$cellTypes == "Endothel"] <- "endothelial cell"
cohen_lung$cellTypes[cohen_lung$cellTypes == "Mon"] <- "monocyte"
cohen_lung$cellTypes[cohen_lung$cellTypes == "Epithel"] <- "epithelial cell"
cohen_lung$cellTypes[cohen_lung$cellTypes %in% c("MacI", "MacII", "MacIII")] <- "macrophage"



tabulaMuris_dropseq_lung$cellTypes <- tabulaMuris_dropseq_lung$cell_ontology_class
tabulaMuris_facs_lung$cellTypes <- tabulaMuris_facs_lung$cell_ontology_class

# exprsMat_mca <- logcounts(mca_lung)

rownames(tabulaMuris_dropseq_lung) <- toupper(rownames(tabulaMuris_dropseq_lung))
rownames(tabulaMuris_facs_lung) <- toupper(rownames(tabulaMuris_facs_lung))



exprsMat_tm <- logcounts(tabulaMuris_dropseq_lung)
exprsMat_tm_facs <- logcounts(tabulaMuris_facs_lung)
exprsMat_cohen <- as.matrix(logcounts(cohen_lung))

rownames(exprsMat_cohen) <- toupper(rownames(exprsMat_cohen))


library(scdney)
library(hopach)
rscript <- list.files("scClassify/")
sapply(rscript, function(x){
  source(paste("scClassify/", x, sep = ""))
})



trainCohenRes <- scClassify(exprsMat_train = exprsMat_cohen,
                            cellTypes_train = cohen_lung$cellTypes,
                            exprsMat_test = list(
                              tm = exprsMat_tm,
                              tm_facs = exprsMat_tm_facs),
                            cellTypes_test = list(
                              cohen = cohen_lung$cellTypes,
                              tm = tabulaMuris_dropseq_lung$cellTypes,
                              tm_facs = tabulaMuris_facs_lung$cellTypes),
                            similarity = c("pearson"),
                            selectFeatures = c("limma"),
                            algorithm = c("WKNN"),
                            parallel = F,
                            ncores = 5,
                            prob_threshold = 0.7,
                            verbose = T,
                            k = 10,
                            topN = 50)



saveRDS(trainCohenRes, file = "results/trainCohenRes_lung.rds")

