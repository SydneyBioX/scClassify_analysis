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
  
  
  
  return(tt)
  
  
}




getEnsembleError <- function(res){
  ensembleErr <- do.call(rbind, lapply(res, function(x) table(x$classifyRes)/length(x$classifyRes)))
  return(ensembleErr)
}


getEnsembleRes <- function(res, exclude = NULL, weight = NULL){
  
  if(!is.null(exclude)){
    keep_method <- !grepl(paste(exclude,collapse = "|"), names(res))
    res <- res[keep_method]
    weight <- weight[keep_method]
  }
  
  if(is.null(weight)){
    weight <- rep(1, length(res))
  }
  
  ensembleResMat <- do.call(cbind, lapply(res, function(x) x$predRes))
  
  ensembleRes <- apply(ensembleResMat, 1, function(x) {
    names(x) <- colnames(ensembleResMat)
    # keep <- x!="unassigned"
    keep <- rep(TRUE, length(x))
    if(sum(keep) == 0){
      data.frame(cellTypes = "unassigned", scores = 0)
    }else{
      getResByWeights(x[keep], weight[keep])
    }
  })
  
  ensembleRes <- do.call(rbind, ensembleRes)
  return(ensembleRes)
}


getResByWeights <- function(res, weight) {
  resType <- unique(res)
  if(length(resType) == 1){
    final <- resType
    scores <- 1
  }else{
    
    
    mat <-sapply(1:length(resType), function(i){
      ifelse(res %in% resType[i], 1, 0) * weight
    })
    colnames(mat) <- resType
    rownames(mat) <- names(res)
    mat_colMeans <- colMeans(mat)
    final <- names(mat_colMeans)[which.max(mat_colMeans)]
    scores <- max(mat_colMeans)
  }
  return(data.frame(cellTypes = final, scores = scores))
}



alpha <- function(e) {
  log((1 - e)/e)
}


ClassifyError2 <- function(cellTypes_pred, cellTypes_test, cellTypes_train){
  if(length(cellTypes_pred)!=length(cellTypes_test)){
    stop("wrong input")
  }
  train_ref <- unique(cellTypes_train)
  test_ref <- unique(cellTypes_test)
  res <- sapply(1:length(cellTypes_pred), function(i){
    if(cellTypes_test[i]%in%train_ref){
      if(cellTypes_pred[i] %in% c("unassigned", "Unassigned")){
        "incorrectly unassigned"
      }else if(cellTypes_pred[i] == "intermediate"){
        "intermediate"
      }else{
        if(cellTypes_test[i] == cellTypes_pred[i]){
          "correct"
        }else if(length(strsplit(cellTypes_pred[i],"_")[[1]])>=2|grepl("Node", cellTypes_pred[i])){
          "intermediate"
        }
        else{
          "misclassified"
        }
      }
    }else{
      if(cellTypes_pred[i] %in% c("unassigned","Unassigned")){
        "correctly unassigned"
      }else{
        "error assigned"
      }
    }
  })
  return(res)
}

theme_yx <- function() {
  theme_bw() +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) +
    theme(plot.title = element_text(hjust = 0.5),
          text = element_text(size=12),
          panel.border = element_rect(colour = "black", fill=NA, size=1.2)) 
  
}

