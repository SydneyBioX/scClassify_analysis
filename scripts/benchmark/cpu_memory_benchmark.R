
.libPaths("/dora/nobackup/yuec/R")



tabulaMuris_facs <- readRDS("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris/00_facs_raw_data/tabulaMuris_facs.rds")
tabulaMuris_dropseq <- readRDS("/albona/nobackup/biostat/datasets/singlecell/tabulaMuris/01_droplet_raw_data/tabulaMuris_dropseq.rds")

tabulaMuris_facs <- tabulaMuris_facs[rowSums(counts(tabulaMuris_facs))!= 0,]
#40000 cells

tabulaMuris_dropseq <- tabulaMuris_dropseq[rowSums(counts(tabulaMuris_dropseq)) != 0,]
#54000 cells

common_hvg <- intersect(rownames(tabulaMuris_dropseq), rownames(tabulaMuris_facs))

exprsMat_cbind <- cbind( tabulaMuris_facs[common_hvg, ],
                         tabulaMuris_dropseq[common_hvg, ])

dim(exprsMat_cbind)

index <- sample(1:ncol(exprsMat_cbind), 48000) 
train_cbind <- exprsMat_cbind[, index   ]
test_cbind <- exprsMat_cbind[,  -c(index)  ] 


sum(table(train_cbind$cell_ontology_class)[order(table(train_cbind$cell_ontology_class), decreasing = T)][1:16])
#35898
sum(table(test_cbind$cell_ontology_class)[order(table(test_cbind$cell_ontology_class), decreasing = T)][1:16])
#36373


top_celltype <- names(table(train_cbind$cell_ontology_class)[order(table(train_cbind$cell_ontology_class), decreasing = T)][1:16])
#select the top 16 cell types




train_cbind   <-  train_cbind[, train_cbind$cell_ontology_class %in% top_celltype]
test_cbind   <-  test_cbind [, test_cbind$cell_ontology_class %in% top_celltype]

dim(train_cbind)
dim(test_cbind)


table(as.character(train_cbind$cell_ontology_class))
table(as.character(test_cbind$cell_ontology_class))





########### set up the benchmark datasets  ######

# needs to be stratified by cell type

set.seed(1)

train_strat <- data.frame(index = 1: ncol(train_cbind), 
                     celltype = train_cbind$cell_ontology_class)


## Taking the sample
library(splitstackshape)
set.seed(1)

strata_train <- function(n){
  out <- stratified(train_strat, c("celltype"), as.integer(n/16) )
  addon <- n - as.integer(n/16)*16
  out = rbind(out, train_strat[1:addon, ])
  
  facs_strata <- train_cbind[, out$index]
  facs_strata_label <- as.character(facs_strata$cell_ontology_class)
  write.csv(logcounts(facs_strata) , paste0("train_", n , ".csv" ))
  write.csv(facs_strata_label , paste0("train_", n , "_label.csv" ))         
  
}

# train - 100 
strata_train(100)

# train - 200
strata_train(200)

# train -  500
strata_train(500)

# train - 1000
strata_train(1000)

# train - 2000
strata_train(2000)

# train - 5000
strata_train(5000)

# train - 10000
strata_train(10000)

# train - 20000
temp <- train_cbind[, sample(1:ncol(train_cbind), 20000, replace = F)]
write.csv(logcounts(temp) , paste0("train_", 20000 , ".csv" ))
write.csv(temp$cell_ontology_class , paste0("train_", 20000 , "_label.csv" ))         

# train - 30000
temp <- train_cbind[, sample(1:ncol(train_cbind), 30000, replace = F)]
write.csv(logcounts(temp) , paste0("train_", 30000 , ".csv" ))
write.csv(temp$cell_ontology_class , paste0("train_", 30000 , "_label.csv" ))         








test_strat <- data.frame(index = 1: ncol(test_cbind), 
                          celltype = test_cbind$cell_ontology_class)

strata_test <- function(n){
  out <- stratified(test_strat, c("celltype"), as.integer(n/16))
  addon <- n - as.integer(n/16)*16
  out = rbind(out, test_strat[1:addon, ])
  
  dropseq_strata <- test_cbind[, out$index]
  dropseq_strata_label <- as.character(dropseq_strata$cell_ontology_class)
  write.csv(logcounts(dropseq_strata) , paste0("test_", n , ".csv" ))
  write.csv(dropseq_strata_label , paste0("test_", n , "_label.csv" ))         
}

# test - 100
strata_test(100)

# test - 200
strata_test(200)

# test - 500
strata_test(500)

# test - 1000
strata_test(1000)

# test - 2000
strata_test(2000)

# test - 5000
strata_test(5000)


# test - 10000
temp <- test_cbind[, sample(1:ncol(test_cbind), 10000, replace = F)]
write.csv(logcounts(temp) , paste0("test_", 10000 , ".csv" ))
write.csv(temp$cell_ontology_class , paste0("test_", 10000 , "_label.csv" ))         

# test - 20000
temp <- test_cbind[, sample(1:ncol(test_cbind), 20000, replace = F)]
write.csv(logcounts(temp) , paste0("test_", 20000 , ".csv" ))
write.csv(temp$cell_ontology_class , paste0("test_", 20000 , "_label.csv" ))         


# test - 30000
temp <- test_cbind[, sample(1:ncol(test_cbind), 30000, replace = F)]
write.csv(logcounts(temp) , paste0("test_", 30000 , ".csv" ))
write.csv(temp$cell_ontology_class , paste0("test_", 30000 , "_label.csv" ))         






 
#########  set up dataset for varying  cell type #########


table(as.character(train_cbind$cell_ontology_class))
# train - 5000, 2 cell type

i_1 <- which(train_cbind$cell_ontology_class  == "B cell")[1:2500]
i_2 <- which(train_cbind$cell_ontology_class  == "basal cell of epidermis")[1:2500]
temp <- train_cbind[, c(i_1, i_2)]
dim(temp)
table(as.character(temp$cell_ontology_class))
write.csv(logcounts(temp) , paste0("train_", 2 , "class.csv" ))
write.csv(temp$cell_ontology_class , paste0("train_", 2 , "class_label.csv" ))         


# train - 5000, 4 cell type
i_1 <- which(train_cbind$cell_ontology_class  == "B cell")[1:2000]
i_2 <- which(train_cbind$cell_ontology_class  == "basal cell of epidermis")[1:1500]
i_3 <- which(train_cbind$cell_ontology_class  == "epithelial cell of large intestine")[1:500]
i_4 <- which(train_cbind$cell_ontology_class  == "hematopoietic stem cell")[1:1000]

temp <- train_cbind[, c(i_1, i_2, i_3, i_4)]
dim(temp)
table(as.character(temp$cell_ontology_class))


write.csv(logcounts(temp) , paste0("train_", 4 , "class.csv" ))
write.csv(temp$cell_ontology_class , paste0("train_", 4 , "class_label.csv" ))         



# train - 5000, 6 cell type
i_1 <- which(train_cbind$cell_ontology_class  == "B cell")[1:1000]
i_2 <- which(train_cbind$cell_ontology_class  == "basal cell of epidermis")[1:1000]
i_3 <- which(train_cbind$cell_ontology_class  == "epithelial cell of large intestine")[1:500]
i_4 <- which(train_cbind$cell_ontology_class  == "hematopoietic stem cell")[1:1000]
i_5 <- which(train_cbind$cell_ontology_class  == "kidney tubule cell")[1:1000]
i_6 <- which(train_cbind$cell_ontology_class  == "macrophage")[1:500]


temp <- train_cbind[, c(i_1, i_2, i_3, i_4, i_5, i_6)]
dim(temp)
table(as.character(temp$cell_ontology_class))

write.csv(logcounts(temp) , paste0("train_", 6 , "class.csv" ))
write.csv(temp$cell_ontology_class , paste0("train_", 6 , "class_label.csv" ))         



# train - 5000, 8 cell type
i_1 <- which(train_cbind$cell_ontology_class  == "B cell")[1:625]
i_2 <- which(train_cbind$cell_ontology_class  == "basal cell of epidermis")[1:625]
i_3 <- which(train_cbind$cell_ontology_class  == "epithelial cell of large intestine")[1:625]
i_4 <- which(train_cbind$cell_ontology_class  == "hematopoietic stem cell")[1:625]
i_5 <- which(train_cbind$cell_ontology_class  == "kidney tubule cell")[1:625]
i_6 <- which(train_cbind$cell_ontology_class  == "macrophage")[1:625]
i_7 <- which(train_cbind$cell_ontology_class  == "basal cell")[1:625]
i_8 <- which(train_cbind$cell_ontology_class  == "endothelial cell")[1:625]


temp <- train_cbind[, c(i_1, i_2, i_3, i_4, i_5, i_6, i_7, i_8)]
dim(temp)
table(as.character(temp$cell_ontology_class))


write.csv(logcounts(temp) , paste0("train_", 8 , "class.csv" ))
write.csv(temp$cell_ontology_class , paste0("train_", 8 , "class_label.csv" ))         






# train - 5000, 10 cell type

i_1 <- which(train_cbind$cell_ontology_class  == "B cell")[1:500]
i_2 <- which(train_cbind$cell_ontology_class  == "basal cell of epidermis")[1:500]
i_3 <- which(train_cbind$cell_ontology_class  == "epithelial cell of large intestine")[1:500]
i_4 <- which(train_cbind$cell_ontology_class  == "hematopoietic stem cell")[1:500]
i_5 <- which(train_cbind$cell_ontology_class  == "kidney tubule cell")[1:500]
i_6 <- which(train_cbind$cell_ontology_class  == "macrophage")[1:500]
i_7 <- which(train_cbind$cell_ontology_class  == "basal cell")[1:500]
i_8 <- which(train_cbind$cell_ontology_class  == "endothelial cell")[1:500]
i_9 <- which(train_cbind$cell_ontology_class  == "fibroblast")[1:500]
i_10 <- which(train_cbind$cell_ontology_class  == "keratinocyte")[1:500]

temp <- train_cbind[, c(i_1, i_2, i_3, i_4, i_5, i_6, i_7, i_8, i_9, i_10)]
dim(temp)
table(as.character(temp$cell_ontology_class))


write.csv(logcounts(temp) , paste0("train_", 10 , "class.csv" ))
write.csv(temp$cell_ontology_class , paste0("train_", 10 , "class_label.csv" ))         




# train - 5000, 12 cell type


i_1 <- which(train_cbind$cell_ontology_class  == "B cell")[1:400]
i_2 <- which(train_cbind$cell_ontology_class  == "basal cell of epidermis")[1:400]
i_3 <- which(train_cbind$cell_ontology_class  == "epithelial cell of large intestine")[1:400]
i_4 <- which(train_cbind$cell_ontology_class  == "hematopoietic stem cell")[1:400]
i_5 <- which(train_cbind$cell_ontology_class  == "kidney tubule cell")[1:400]
i_6 <- which(train_cbind$cell_ontology_class  == "macrophage")[1:400]
i_7 <- which(train_cbind$cell_ontology_class  == "basal cell")[1:400]
i_8 <- which(train_cbind$cell_ontology_class  == "endothelial cell")[1:400]
i_9 <- which(train_cbind$cell_ontology_class  == "fibroblast")[1:400]
i_10 <- which(train_cbind$cell_ontology_class  == "keratinocyte")[1:400]
i_11 <- which(train_cbind$cell_ontology_class  == "leukocyte")[1:400]
i_12 <- which(train_cbind$cell_ontology_class  == "mesenchymal cell")[1:600]


temp <- train_cbind[, c(i_1, i_2, i_3, i_4, i_5, i_6, i_7, i_8, i_9, i_10 , i_11 , i_12)]
dim(temp)
table(as.character(temp$cell_ontology_class))


write.csv(logcounts(temp) , paste0("train_", 12 , "class.csv" ))
write.csv(temp$cell_ontology_class , paste0("train_", 12 , "class_label.csv" ))         













#### to measure performance #############

#tools like 
#Rprof,
#Rprofmem, 
#profmem::profmem,
#bench::mark or 



smallest.sv <- function(){
  A <- matrix(rnorm(1e6), 1e1);
  mysvd <- svd(A);
 
  return(tail(mysvd$d, 1));
}





#Using Rprofmem (Enable momory profiling is a compile-time option: ./configure --enable-memory-profiling)
#best method for parallel code 
Rprofmem("Rprofmem.out")
x <- smallest.sv()
Rprofmem(NULL)
sum(as.numeric(read.table("Rprofmem.out", comment.char = ":", sep= ",")[,1]), na.rm=TRUE)


# best method
#88101752
#Writes out them memory amount when it is allocated
library(profmem) 
total(profmem(smallest.sv()))



# best method, (wrapper around profmem)  but does not work for parallel code 
library(bench) 
mark(smallest.sv())[, "mem_alloc"]
mark(smallest.sv() , time_unit = "s")[, "total_time"]
#84MB
#Warning message:
#Some expressions had a GC in every iteration; so filtering is disabled. 

#useless
library(profvis)
profvis(x <- smallest.sv())
#opens a browser window where you can read under Memory -23.0 | 45.9


