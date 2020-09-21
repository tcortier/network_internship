require(fields)
require(parallel)
require(stringr)

check_boundaries_2 <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x <     zero] <-     zero
  x
}

LBMbar <- function(X) {
  X.bar <- 1 - X 
  return(X.bar)
}

.softmax <- function(x) {
  b <- max(x)
  exp(x - b) / sum(exp(x - b))
}

check_boundaries <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x > 1 - zero] <- 1 - zero
  x[x <     zero] <-     zero
  return(x)
}

quad_form <- function(A,x) {t(x) %*% A %*% x}

logistic <- function(x) {1/(1 + exp(-x))}
logit    <- function(x) {log(x/(1 - x))}

clustering_indicator <- function(clustering) {##Returns the tau variable from a vector containing the clustering data
  nBlocks <- length(unique(clustering))
  nNodes  <- length(clustering)
  Z <- matrix(0,nNodes, nBlocks)
  Z[cbind(seq.int(nNodes), clustering)] <- 1
  return(Z)
}



clustinit_LBM<-function(connectivity,Q1,Q2,type="hierarchical_clust"){
  ### Returns a list of two vectors containing which cluster belongs to which node from a defined number of nodes
  ### The three possible types are "hierarchical_clust", "spectral_clust","kmeans_clust"
  if (type=="hierarchical_clust"){
    if (Q1 > 1) {
      D <- as.matrix(dist(connectivity, method = "manhattan"))
      D <- as.dist(ape::additive(D))
      cl01 <- cutree(hclust(D, method = "ward.D2"), Q1)
    } else {
      cl01 <- rep(1L,nrow(connectivity))
    }
    if (Q2 > 1) {
      D <- as.matrix(dist(t(connectivity), method = "manhattan"))
      D <- as.dist(ape::additive(D))
      cl02 <- cutree(hclust(D, method = "ward.D2"), Q2)
    } else {
      cl02 <- rep(1L,nrow(t(connectivity)))
    }
    cl0=list(cl01,cl02)
    return(cl0)
  }
  else if (type=="spectral_clust"){
    n1=dim(connectivity)[1]
    n2=dim(connectivity)[2]
    if (Q1 > 1) {
      degrees = apply(connectivity, 1,sum)
      degrees = (1/degrees)%x%matrix(1,1,n2)
      L=degrees*connectivity
      dec=svd(L)
      U = dec$u[,1:Q1]
      km=kmeans(U,Q1)
      tau=rdist(km$centers,U)
      cl01=apply(tau,2,which.max)
    }
    else{
      cl01 <- rep(1L,nrow(connectivity))
    }
    if (Q2 > 1) {
      degrees = apply(connectivity, 2,sum)
      degrees = (1/degrees)%x%matrix(1,1,n1)
      L=degrees*t(connectivity)
      dec=svd(L)
      U = dec$u[,1:Q2]
      km=kmeans(U,Q2)
      tau=rdist(km$centers,U)
      cl02=apply(tau,2,which.max)
    }
    else{
      cl02 <- rep(1L,ncol(connectivity))
    }
    cl0=list(cl01,cl02)
    return(cl0)
  }
  else if (type=="kmeans_clust"){
    if (Q1 > 1) {
      D  <- as.matrix(dist(connectivity, method = "euclidean"))
      cl01 <- as.integer(kmeans(ape::additive(D), Q1, nstart = 50, iter.max = 100)$cl)
    } else {
      cl01 <- rep(1L, nrow(connectivity))
    }
    if (Q2 > 1) {
      D  <- as.matrix(dist(t(connectivity), method = "euclidean"))
      cl02 <- as.integer(kmeans(ape::additive(D), Q2, nstart = 50, iter.max = 100)$cl)
    } else {
      cl02 <- rep(1L, nrow(t(connectivity)))
    }
    cl0=list(cl01,cl02)
    return(cl0)
  }
}

membertoclust<-function(membership){
  return(apply(membership,1,which.max))
}

LBM_plot<-function(models){
  ##Function that plots ICLs according to the number of groups 
  modnames=names(models)
  ICL=unlist(lapply(models, '[[','ICL'))
  Q=as.numeric(str_sub(str_extract(modnames, ".&|..&"),1,-2))+as.numeric(str_sub(str_extract(modnames, "&.|&.."),2,-1))
  plot(Q,ICL,xlab="Q1+Q2",ylab="ICL")
}

forward_explo<-function(models,k1,k2,connectivity,param){
  ##Function to explore all possible sharing combinations of groups from k1 and k2 groups
  cl01<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="&")]]$membership1)
  cl02<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="&")]]$membership2)
  n1=length(unique(cl01))
  n2=length(unique(cl02))
  print(n1)
  print(n2)
  best_one=list()
  if (n1==k1&n2==k2){
    candidates <- mclapply(1:n1, function(j) {
      cl1 <- cl01
      J  <- which(cl1 == j)
      if (length(J) > 1) {
        J1 <- base::sample(J, floor(length(J)/2))
        J2 <- setdiff(J, J1)
        cl1[J1] <- j; 
        cl1[J2] <- n1 + 1
        model=LBM_VEM(connectivity,n1+1,k2,clustering_indicator(cl1),clustering_indicator(cl02),param)
      }
      else {
        model=models[[paste(as.character(k1+1),as.character(k2),sep="&")]]
      }
      return(model)
    },mc.cores = param$cores)
    best_one[[1]] <- candidates[[which.min(sapply(candidates, function(candidate) candidate$ICL))]]
    
    candidates <- mclapply(1:n2, function(j) {
      cl2 <- cl02
      J  <- which(cl2 == j)
      if (length(J) > 1) {
        J1 <- base::sample(J, floor(length(J)/2))
        J2 <- setdiff(J, J1)
        cl2[J1] <- j; 
        cl2[J2] <- n2 + 1
        model=LBM_VEM(connectivity,k1,k2+1,clustering_indicator(cl01),clustering_indicator(cl2),param)
      }
      else {
        model=models[[paste(as.character(k1),as.character(k2+1),sep="&")]]
      }
      return(model)
    },mc.cores = param$cores)
    best_one[[2]] <- candidates[[which.min(sapply(candidates, function(candidate) candidate$ICL))]]
  }
  else{
    best_one[[1]]=models[[paste(as.character(k1+1),as.character(k2),sep="&")]]
    best_one[[2]]=models[[paste(as.character(k1),as.character(k2+1),sep="&")]]
  }
  return(best_one)
}



backward_explo<-function(models,k1,k2,connectivity,param){
  cl01<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="&")]]$membership1)
  cl02<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="&")]]$membership2)
  cl1 <- factor(cl01)
  cl2 <- factor(cl02)
  n1=nlevels(cl1)
  n2=nlevels(cl2)
  print(n1)
  print(n2)
  best_one=list()
  if (n1==k1&n2==k2){
    candidates<-mclapply(combn(n1, 2, simplify = FALSE),function(couple){
      cl_fusion1 <- cl1
      levels(cl_fusion1)[which(levels(cl_fusion1) == paste(couple[1]))] <- paste(couple[2])
      levels(cl_fusion1) <- as.character(1:(n1-1))
      cl_fusion1<-as.numeric(cl_fusion1)
      model=LBM_VEM(connectivity,n1-1,k2,clustering_indicator(cl_fusion1),clustering_indicator(cl02),param)
      return(model)
    },mc.cores = param$cores)
    best_one[[1]] <- candidates[[which.min(sapply(candidates, function(candidate) candidate$ICL))]]
    
    candidates<-mclapply(combn(n2, 2, simplify = FALSE),function(couple){
      cl_fusion2 <- cl2
      levels(cl_fusion2)[which(levels(cl_fusion2) == paste(couple[1]))] <- paste(couple[2])
      levels(cl_fusion2) <- as.character(1:(n2-1))
      cl_fusion2<-as.numeric(cl_fusion2)
      model=LBM_VEM(connectivity,k1,n2-1,clustering_indicator(cl01),clustering_indicator(cl_fusion2),param)
      return(model)
    },mc.cores = param$cores)
    best_one[[2]] <- candidates[[which.min(sapply(candidates, function(candidate) candidate$ICL))]]
  }
  else{
    best_one[[1]]=models[[paste(as.character(k1-1),as.character(k2),sep="&")]]
    best_one[[2]]=models[[paste(as.character(k1),as.character(k2-1),sep="&")]]
  }
  return(best_one)
}



LBM_update_tau<-function(connectivity,alpha1,alpha2,pi,tau1,tau2){##Update of tau1 and tau2 taking into account sampling for a binomiale distribution
  barconnectivity=LBMbar(connectivity)
  Q1=length(alpha1)
  Q2=length(alpha2)
  if ((Q1>1)&(Q2>1)){
    newtau1<-connectivity %*% tau2 %*% t(log(pi)) + barconnectivity %*% tau2 %*% t(log(1 - pi))
    newtau1 <- t(apply(sweep(newtau1, 2, log(alpha1), "+"), 1, .softmax))
    newtau2<-t(connectivity) %*% tau1 %*% log(pi) + t(barconnectivity) %*% tau1 %*% log(1 - pi)
    newtau2 <- t(apply(sweep(newtau2, 2, log(alpha2), "+"), 1, .softmax))
  }
  else if ((Q1>1)&(Q2==1)){
    newtau1<-connectivity %*% tau2 %*% t(log(pi)) + barconnectivity %*% tau2 %*% t(log(1 - pi))
    newtau1 <- t(apply(sweep(newtau1, 2, log(alpha1), "+"), 1, .softmax))
    newtau2=tau2
  }
  else if ((Q1==1)&(Q2>1)){
    newtau2<-t(connectivity) %*% tau1 %*% log(pi) + t(barconnectivity) %*% tau1 %*% log(1 - pi)
    newtau2 <- t(apply(sweep(newtau2, 2, log(alpha2), "+"), 1, .softmax))
    newtau1=tau1
  }
  else if ((Q1==1)&(Q2==1)){
    newtau1=tau1
    newtau2=tau2
  }
  return(list(newtau1,newtau2))
}

LBM_update_tau_poisson<-function(connectivity,alpha1,alpha2,lambda,tau1,tau2){##Update of tau1 and tau2 taking into account sampling for a poisson distribution
  Q1=length(alpha1)
  Q2=length(alpha2)
  N1=dim(connectivity)[1]
  N2=dim(connectivity)[2]
  if ((Q1>1)&(Q2>1)){
    newtau1<-connectivity %*% tau2 %*% t(log(lambda)) - matrix(1,N1,N2) %*% tau2 %*% t(lambda)-log(factorial(connectivity))%*%tau2%*%matrix(1,Q2,Q1)
    newtau1 <- t(apply(sweep(newtau1, 2, log(alpha1), "+"), 1, .softmax))
    newtau2<-t(connectivity) %*% tau1 %*% log(lambda) - matrix(1,N2,N1)%*% tau1 %*% lambda - t(log(factorial(connectivity)))%*%tau1%*%matrix(1,Q1,Q2)
    newtau2 <- t(apply(sweep(newtau2, 2, log(alpha2), "+"), 1, .softmax))
  }
  else if ((Q1>1)&(Q2==1)){
    newtau1<-connectivity %*% tau2 %*% t(log(lambda)) - matrix(1,N1,N2) %*% tau2 %*% t(lambda)-log(factorial(connectivity))%*%tau2%*%matrix(1,Q2,Q1)
    newtau1 <- t(apply(sweep(newtau1, 2, log(alpha1), "+"), 1, .softmax))
    newtau2=tau2
  }
  else if ((Q1==1)&(Q2>1)){
    newtau2<-t(connectivity) %*% tau1 %*% log(lambda) - matrix(1,N2,N1)%*% tau1 %*% lambda - t(log(factorial(connectivity)))%*%tau1%*%matrix(1,Q1,Q2)
    newtau2 <- t(apply(sweep(newtau2, 2, log(alpha2), "+"), 1, .softmax))
    newtau1=tau1
  }
  else if ((Q1==1)&(Q2==1)){
    newtau1=tau1
    newtau2=tau2
  }
  return(list(newtau1,newtau2))
}


LBM_update_pi = function(connectivity,tau1,tau2) { ##We update pi for a binomial distribution
  N1<-dim(connectivity)[1]
  N2<-dim(connectivity)[2]
  pi    <- check_boundaries(t(tau1)%*%connectivity%*%tau2 / t(tau1)%*%matrix(1,N1,N2)%*%tau2)
  return(pi)
}

LBM_update_lambda = function(connectivity,tau1,tau2) { ##Update lambda for a poisson distribution
  N1<-dim(connectivity)[1]
  N2<-dim(connectivity)[2]
  lambda    <- check_boundaries_2(t(tau1)%*%connectivity%*%tau2 / t(tau1)%*%matrix(1,N1,N2)%*%tau2)
  return(lambda)
}

LBM_update_alpha = function(tau) { ##Update alpha
  alpha <- check_boundaries(colMeans(tau))
  return(alpha)
}

memberships = function(tau) {##Construct Z estimate (0 or 1) from tau
  dimtau=dim(tau)
  mem=matrix(0,dimtau[1],dimtau[2])
  for (i in 1:dimtau[1]){
    a=which.max(tau[i,])
    mem[i,a]=1
  }
  return(mem)
}

LBM_ICL<-function(members1,members2,alpha1,alpha2,pi, connectivity){##Compute ICL from model's data for bernoulli distribution
  Q1=length(alpha1)
  Q2=length(alpha2)
  N1=dim(connectivity)[1]
  N2=dim(connectivity)[2]
  barconnectivity=LBMbar(connectivity)
  logL=sum(members1%*%log(alpha1))+sum(members2%*%log(alpha2))+sum((t(members1)%*%connectivity%*%members2)*log(pi))+sum((t(members1)%*%barconnectivity%*%members2)*log(1-pi))
  pena=log(N1)*(Q1-1)/2+log(N1*(N1+1)/2)*Q1*(Q1+1)/4+log(N2)*(Q2-1)/2+log(N2*(N2+1)/2)*Q2*(Q2+1)/4
  return(logL-pena)
}

LBM_ICL_poisson<-function(members1,members2,alpha1,alpha2,lambda, connectivity){##Compute ICL from model's data for poisson distribution
  Q1=length(alpha1)
  Q2=length(alpha2)
  N1=dim(connectivity)[1]
  N2=dim(connectivity)[2]
  logL=sum(members1%*%log(alpha1))+sum(members2%*%log(alpha2))+sum((t(members1)%*%connectivity%*%members2)*log(lambda))-sum((t(members1)%*%matrix(1,N1,N2)%*%members2)*lambda)-sum((t(members1)%*%log(factorial(connectivity))%*%members2))
  pena=log(N1)*(Q1-1)/2+log(N1*(N1+1)/2)*Q1*(Q1+1)/4+log(N2)*(Q2-1)/2+log(N2*(N2+1)/2)*Q2*(Q2+1)/4
  return(logL-pena)
}


LBM_VEM<-function(connectivity,Q1,Q2,tau1=c(),tau2=c(),param) {
  ##param=list(type=Bernoulli,maxIter=50,fixPointIter=3,threshold=1e-3,trace=TRUE,cores=1)
  ##type in param reprensents the distribution to use in the model
  ##maxIter is the maximum number of iteraction before the algorithlm stop even if the convergence conditions are not reached
  ##fixPointIter is the number of iteraction to use in the fixed point algorithm for tau
  ##threshold is The convergence threshold
  ##trace is the parameter to display some data of the algorithm
  ##cores how many cores to be used
  if (param$type=='bernoulli'){
    if (param$trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
    ## Initialization of quantities that monitor convergence
    delta     <- vector("numeric", param$maxIter)
    pi=matrix(0.5,Q1,Q2)
    i <- 0; 
    cond <- FALSE
    
    while (!cond) {
      
      i <- i + 1
      
      ## ______________________________________________________
      ## M-step
      #
      # update the parameters of the SBM (a.k.a alpha and pi)
      alpha1=LBM_update_alpha(tau1)
      alpha2=LBM_update_alpha(tau2)
      new_pi=LBM_update_pi(connectivity,tau1,tau2)
      
      ## ______________________________________________________
      ## Variational E-Step
      #
      for (k in seq.int(param$fixPointIter)) {
        tau=LBM_update_tau(connectivity,alpha1,alpha2,new_pi,tau1,tau2)
        tau1=tau[[1]]
        tau2=tau[[2]]
      }
      
      ## Check convergence
      delta[i] <- sqrt(sum((new_pi - pi)^2)) / sqrt(sum((pi)^2))
      cond     <- (i > param$maxIter) |  (delta[i] < param$threshold)
      pi=new_pi
    }
    members1=memberships(tau1)
    members2=memberships(tau2)
    icl=LBM_ICL(members1,members2,alpha1,alpha2,pi,connectivity)
    res=list()
    res$tau1=tau1
    res$tau2=tau2
    res$pi=pi
    res$alpha1=alpha1
    res$alpha2=alpha2
    res$ICL=icl
    res$membership1=members1
    res$membership2=members2
    return(res)
  }
  else if(param$type=='poisson'){
    if (param$trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
    ## Initialization of quantities that monitor convergence
    delta     <- vector("numeric", param$maxIter)
    lambda=matrix(0,Q1,Q2)
    i <- 0; 
    cond <- FALSE
    
    while (!cond) {
      
      i <- i + 1
      
      ## ______________________________________________________
      ## M-step
      #
      # update the parameters of the SBM (a.k.a alpha and pi)
      alpha1=LBM_update_alpha(tau1)
      alpha2=LBM_update_alpha(tau2)
      new_lambda=LBM_update_lambda(connectivity,tau1,tau2)
      
      ## ______________________________________________________
      ## Variational E-Step
      #
      for (k in seq.int(param$fixPointIter)) {
        tau=LBM_update_tau_poisson(connectivity,alpha1,alpha2,new_lambda,tau1,tau2)
        tau1=tau[[1]]
        tau2=tau[[2]]
      }
      ## Check convergence
      delta[i] <- sqrt(sum((new_lambda - lambda)^2)) / sqrt(sum((lambda)^2))
      cond     <- (i > param$maxIter) |  (delta[i] < param$threshold)
      lambda=new_lambda
    }
    members1=memberships(tau1)
    members2=memberships(tau2)
    icl=LBM_ICL_poisson(members1,members2,alpha1,alpha2,lambda,connectivity)
    res=list()
    res$tau1=tau1
    res$tau2=tau2
    res$lambda=lambda
    res$alpha1=alpha1
    res$alpha2=alpha2
    res$ICL=icl
    res$membership1=members1
    res$membership2=members2
    return(res)
  }
}

