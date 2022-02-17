#rmtl_X: list of feature matrices, each matrix is nk*p
#rmtl_Y: list of gene expression, each matrix is nk*1
#type: Sparse (MTL with L1) Sparse_L21 (MTL with L1 and L21)
SMTL <- function(rmtl_X,rmtl_Y,type = "SParse_L21"){
  if(type == "SParse_L21"){
    K = length(rmtl_X)
    cvfit<-cvSMTL_L21(rmtl_X, rmtl_Y)
    m=admm.iters(Y = rmtl_Y,X = rmtl_X,lambda1 = cvfit$Lam1.min,lambda2   = cvfit$Lam2.min)$theta
    m = matrix(unlist(m),ncol = K,byrow = T)
    
  }else if(type == "Sparse"){
    cvfit<-cvMTL(X = rmtl_X,Y=rmtl_Y,type = "Regression",Regularization = "Lasso")
    m=t(MTL(X = rmtl_X,Y=rmtl_Y,Lam1 = cvfit$Lam1.min,Lam1_seq = cvfit$Lam1_seq)$W)
  }else if(type == "LASSO"){
    K = length(rmtl_X)
    m = t(sapply(1:K,function(i){
      fm <- cv.glmnet(rmtl_X[[i]],y =as.numeric(rmtl_Y[[i]]),intercept = F)
      fm <- glmnet(x=rmtl_X[[i]],y =as.numeric(rmtl_Y[[i]]),lambda = fm$lambda.1se,intercept = F)
      return(coef(fm)[-1,1])
    }))
  }
  colnames(m) = colnames(rmtl_X[[1]])
  #rownames(m) = rownames(rmtl_Y[[1]])
  return(m)
}

admm.iters = function(Y,X,lambda1=.15,lambda2=.2,rho=1,rho.increment=1,weights,maxiter = 1000,tol=1e-5,truncate=1e-5,scale = F){
  K = length(Y)
  for(k in 1:K){Y[[k]]=as.matrix(Y[[k]])}
  P = dim(X[[1]])[2]
  if(scale){
    for(k in 1:K){

      Y[[k]]= scale(Y[[k]])
    }
  }

  
  # initialize theta:
  theta = list(); for(k in 1:K){theta[[k]] = rep(0,P)}
  #for(k in 1:K){theta[[k]] = 1/S[[k]]}
  # initialize Z:
  Q = list(); for(k in 1:K){Q[[k]] = rep(0,P)}
  # initialize U:
  U = list();	for(k in 1:K){U[[k]] = rep(0,P)}
  
  # # initialize lambdas:
  # lam1 = lambda1
  # lam2 = lambda2
  # 
  # iterations:
  iter=0
  diff_value = 10
  while((iter==0) || (iter<maxiter && diff_value > tol))
  {
    # update theta to minimize loss + <S,theta> + rho/2||theta - Q + U ||^2_2:
    theta.prev = theta
    for(k in 1:K)
    {
      #update theta for squared loss
      theta[[k]] = updateTheta(X= X[[k]],Y= Y[[k]],P = diag(P),rho = as.double(rho),U= U[[k]],Q = Q[[k]])
    }
    
    # update Q:
    # define A matrices:
    A = list()
    for(k in 1:K){ A[[k]] = theta[[k]] + U[[k]] }
    
    # minimize rho/2 ||Q-A||_F^2 + P(Q):
    Q = dsgl(A,rho,lambda1,lambda2)
    
    
    # update the dual variable U:
    for(k in 1:K){U[[k]] = U[[k]] + (theta[[k]]-Q[[k]])}
    
    # bookkeeping:
    iter = iter+1
    diff_value = 0
    for(k in 1:K) {diff_value = diff_value + sum(abs(theta[[k]] - theta.prev[[k]])) / sum(abs(theta.prev[[k]]))}
    # increment rho by a constant factor:
    rho = rho * rho.increment
  }
  diff = 0; for(k in 1:K){diff = diff + sum(abs(theta[[k]]-Q[[k]]))}
  
  # round very small theta entries down to zero:
  if(dim(theta[[1]])[1]>0)
  {
    for(k in 1:K)
    {
      rounddown = abs(theta[[k]])<truncate; diag(rounddown)=FALSE
      theta[[k]]=theta[[k]]*(1-rounddown)
    }}
  
  out = list(theta=theta,Q=Q,diff=diff,iters=iter)
  return(out)
}

soft <-function(a,lam){ # if last argument is FALSE, soft-threshold a matrix but don't penalize the diagonal
  out <- sign(a)*pmax(0, abs(a)-lam)
  return(out)
}

dsgl <-function(A,L,lam1,lam2){
  lam1 = lam1*1/L
  lam2 = lam2*1/L
  
  if(is.matrix(A[[1]])) {p=dim(A[[1]])[1]}
  if(is.vector(A[[1]])) {p=length(A[[1]])}
  K=length(A)
  softA = A
  for(k in 1:K) {softA[[k]] = soft(A[[k]],lam1) } 
  normsoftA = A[[1]]*0
  for(k in 1:K) {normsoftA = normsoftA + (softA[[k]])^2}
  
  normsoftA = sqrt(normsoftA)
  
  notshrunk = (normsoftA>lam2)*1
  # reset 0 elements of normsoftA to 1 so we don't get NAs later. 
  normsoftA = normsoftA + (1-notshrunk)
  
  out = A
  for(k in 1:K)
  {
    out[[k]] = softA[[k]]*(1-lam2/normsoftA)
    out[[k]] = out[[k]]*notshrunk
  }
  return(out)
}


