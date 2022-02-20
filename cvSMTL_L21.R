cvSMTL_L21 <- function(X, Y, Lam1_seq=10^seq(1,-4, -1), Lam2_seq=10^seq(2,-4, -1),nfolds=5, ncores=2, parallel=FALSE){
  #test vilidity of input data
  if (!missing(X) & !missing(Y)){
    if (all(sapply(X, class)!="matrix")){
      X <- lapply(X, function(x){as.matrix(x)})
    }
    if (all(sapply(Y, class)!="matrix")){
      Y <- lapply(Y, function(x){as.matrix(x)})
    }
  }else{
    stop("data X or Y doesnot exists")
  }
  task_num <- length(X)
  # if(stratify & type=="Regression"){
  #   stop("stratified CV is not applicable to regression")}
  cvPar <- getCVPartition(Y, nfolds)
  
  #cv
  if (!parallel){
    cvm <- matrix(0, length(Lam1_seq),length(Lam2_seq))
    for (i in 1:nfolds){
      cv_Xtr <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[1]][[x]], ])
      cv_Ytr <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[1]][[x]]])
      cv_Xte <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[2]][[x]], ])
      cv_Yte <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[2]][[x]]])
      #opt <- opts
      for (lam1_idx in 1: length(Lam1_seq)){
        for (lam2_idx in 1: length(Lam2_seq)) {
          m <- admm.iters(Y = cv_Ytr,X=cv_Xtr,lambda1=Lam1_seq[lam1_idx],lambda2=Lam2_seq[lam2_idx])
          cv_err <- calcError(m, newX=cv_Xte, newY=cv_Yte)
          cvm[lam1_idx,lam2_idx] = cvm[lam1_idx,lam2_idx]+cv_err
        }
        # m <- MTL(X=cv_Xtr, Y=cv_Ytr, type=type,
        #          Regularization=Regularization, Lam1=Lam1_seq[lam1_idx],
        #          Lam2=Lam2, opts=opt, k=k, G=G)
        #non sparse model training
        # if (!is.element(Regularization, c("Graph", "CMTL"))){
        #   opt$init <- 1
        #   opt$W0 <- m$W
        #   opt$C0 <- m$C
        # }
        
      }
    }
    cvm = cvm/nfolds
  } else {
    requireNamespace('doParallel')
    requireNamespace('foreach')
    doParallel::registerDoParallel(ncores)
    cvm <- foreach::foreach(i = 1:nfolds, .combine="cbind") %dopar%{
      cv_Xtr <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[1]][[x]], ])
      cv_Ytr <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[1]][[x]]])
      cv_Xte <- lapply(c(1:task_num),
                       function(x) X[[x]][cvPar[[i]][[2]][[x]], ])
      cv_Yte <- lapply(c(1:task_num),
                       function(x) Y[[x]][cvPar[[i]][[2]][[x]]])
      opt <- opts
      cvVec=rep(0, length(Lam1_seq))
      for (lam1_idx in 1: length(Lam1_seq)){
        m <- MTL(X=cv_Xtr, Y=cv_Ytr, type=type,
                 Regularization=Regularization, Lam1=Lam1_seq[lam1_idx],
                 Lam2=Lam2, opts=opt, k=k, G=G)
        #non sparse model training
        if (!is.element(Regularization, c("Graph", "CMTL")) ){
          opt$init <- 1
          opt$W0 <- m$W
          opt$C0 <- m$C
        }
        cv_err <- calcError(m, newX=cv_Xte, newY=cv_Yte)
        cvVec[lam1_idx] <- cv_err
      }
      return(cvVec)
    }
    cvm <- rowMeans(cvm)
  }
  
  best_idx_l1 <- which(cvm == min(cvm), arr.ind = TRUE)[1]
  best_idx_l2 <- which(cvm == min(cvm), arr.ind = TRUE)[2]
  cv <- list(Lam1_seq=Lam1_seq, Lam1.min=Lam1_seq[best_idx_l1],
             Lam2_seq=Lam2_seq, Lam2.min=Lam2_seq[best_idx_l2], cvm=cvm)
  class(cv) <- "cvMTL"
  return(cv)
}



getCVPartition <- function(Y, cv_fold){
  task_num = length(Y);
  
  randIdx <- lapply(Y, function(x) sample(1:length(x),
                                          length(x), replace = FALSE))        
  cvPar = {};
  for (cv_idx in 1: cv_fold){
    # buid cross validation data splittings for each task.
    cvTrain = {};
    cvTest = {};
    
    #stratified cross validation
    for (t in 1: task_num){
      task_sample_size <- length(Y[[t]]);
      
      te_idx <- seq(cv_idx, task_sample_size, by=cv_fold)
      tr_idx <- seq(1,task_sample_size)[!is.element(1:task_sample_size, te_idx)];
      
      cvTrain[[t]] = randIdx[[t]][tr_idx]
      cvTest[[t]] = randIdx[[t]][te_idx]
    }
    
    cvPar[[cv_idx]]=list(cvTrain, cvTest);
  }
  return(cvPar)
}

calcError <- function(m, newX=NULL, newY=NULL){
  # if(class(m)!="MTL"){
  #   stop("The first arguement is not a MTL model")}
  if(!is.null(newX) & !is.null(newY)){
    task_num <- length(newY)
    yhat <- predict.MTL(m,newX)
    error <- sapply(1:task_num, function(x)
      mean((newY[[x]]-yhat[[x]])^2))
    return(mean(error))
  }else{stop(" no new data (X or Y) are provided ")}
}

predict.MTL <- function(object, newX=NULL, ...){
  if(!is.null(newX)){
    task_num <- length(newX)
    score <- lapply(c(1:task_num), function(x)
      newX[[x]] %*% object$theta[[x]] )
    return(score)
  }else{stop("no new data (X) is provided")}
}