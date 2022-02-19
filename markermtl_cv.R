# len: number of genes to be used
marker_smtl_cv = function(ncore = 40,X = rmtl_X,data = stab_umi_lcpm2_sub2,len = 3,parallel= T,type = "SParse_L21") { 
  Y_nam= names(data)
  X = X[Y_nam]
  if(parallel){
    
  
  # , log = T
  cl <- makeCluster(ncore)
  registerDoSNOW(cl)
  
  sig_mat=foreach(i=1:len, .packages = c("MASS","RMTL","psych","corpcor","fields","glmnet"), .errorhandling='pass') %dopar% {
    #result <- tryCatch({
    rmtl_Y = lapply(data, function(x){
      y = as.matrix(x[i,])
      return(y)
    })
    Rcpp::sourceCpp("D:/Manqi/sc-meth/SMTL/SMTL_code/upatetheta.cpp")
    source("D:/Manqi/sc-meth/SMTL/SMTL_code/cvSMTL_L21.R")
    source("D:/Manqi/sc-meth/SMTL/SMTL_code/SMTL_L21.r")
    #m = SMTL(rmtl_X = X, rmtl_Y = rmtl_Y,type = type)
    #cvfit<-cvMTL(X = X,Y=rmtl_Y,type = "Regression",Regularization = "Lasso")
    #m=MTL(X = X,Y=rmtl_Y,Lam1 = cvfit$Lam1.min,Lam1_seq = cvfit$Lam1_seq)
    cvm=cvSMTL_L21(X, rmtl_Y)$cvm
    return(cvm)
   # },
      #  error = function(err) {
      #     return('an error')}) # END tryCatch
    
   # return(result)  # return your result to outputlist
  }
  
  stopCluster(cl)
  }else{
    sig_mat = list()
    for (i in 1:len) {
      rmtl_Y = lapply(data, function(x){
        y = as.matrix(x[i,])
        return(y)
      })

      m = SMTL(rmtl_X = X, rmtl_Y = rmtl_Y,type = type)
      sig_mat[[i]] = m
      
    }
  }
  return(sig_mat)
}

SMTL <- function(rmtl_X,rmtl_Y,type = "SParse_L21"){
  if(type == "SParse_L21"){
    K = length(rmtl_X)
    cvfit<-cvSMTL_L21(rmtl_X, rmtl_Y)
    m=admm.iters(Y = rmtl_Y,X = rmtl_X,lambda1 = cvfit$Lam1.min,lambda2   = cvfit$Lam2.min)$theta
    m = matrix(unlist(m),ncol = K)
  }else if(type == "Sparse"){
    cvfit<-cvMTL(X = rmtl_X,Y=rmtl_Y,type = "Regression",Regularization = "Lasso")
    m=MTL(X = rmtl_X,Y=rmtl_Y,Lam1 = cvfit$Lam1.min,Lam1_seq = cvfit$Lam1_seq)$W
  }
  return(m)
}


# marker_rmtl = function(ncore = 30,X = rmtl_X,data = stab_umi_lcpm2_sub2,len) { 
#   # , log = T
#   #registerDoMC(cores=ncore)
#   cl <- makeCluster(ncore)
#   registerDoSNOW(cl)
#   sig_mat=foreach(i=1:len, .packages = c("MASS","RMTL","psych","corpcor","fields"), .errorhandling='pass') %dopar% {
#     result <- tryCatch({
#       source("D:/Manqi/sc-meth/JGL-master/R/cvSMTL_L21.R")
#       source("D:/Manqi/sc-meth/JGL-master/R/SMTL_L21.r")
#       rmtl_Y = list(darmanis = as.matrix(data$Darmanis[i,]), Habib =as.matrix(data$Habib[i,]),Welch = as.matrix(data$Welch[i,]))
#       cvfit<-cvSMTL_L21(X, rmtl_Y)
#       m=admm.iters(Y = rmtl_Y,X,lambda1 = cvfit$Lam1.min,lambda2   = cvfit$Lam2.min)
#       return(m$theta)
#     },
#     error = function(err) {
#       return('an error')}) # END tryCatch
#     
#     return(result)  # return your result to outputlist
#   }
#   
#   stopCluster(cl)
#   return(sig_mat)
# }