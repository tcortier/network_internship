require(fields)
require(blockmodels)
require(parallel)
require(stats)
require(stringr)

LBMbar <- function(X) {
  X.bar <- 1 - X 
  return(X.bar)
}

bar_nu <- function(connectivity,nu) {
  nubar <- 1 - nu - connectivity
  return(nubar)
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

.logit<-function(x){
  return(log(x/(1-x)))
}

.logistic<-function(x){
  return(1/(1+exp(-x)))
}

quad_form <- function(A,x) {t(x) %*% A %*% x}


clustering_indicator <- function(clustering) {##Returns the tau variable from a vector containing the clustering data
  nBlocks <- length(unique(clustering))
  nNodes  <- length(clustering)
  Z <- matrix(0,nNodes, nBlocks)
  Z[cbind(seq.int(nNodes), clustering)] <- 1
  return(Z)
}


##Essai pour l'initialisation

##Next question
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
      cl02 <- rep(1L,nrow(connectivity))
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
  ##Returns a vector with the cluster to which the nodes belong from tau
  return(apply(membership,1,which.max))
}

LBM_plot<-function(models){
  ##Function that plots ICLs according to the number of groups 
  modnames=names(models)
  ICL=unlist(lapply(models, '[[','ICL'))
  Q=as.numeric(str_sub(str_extract(modnames, ".&|..&"),1,-2))+as.numeric(str_sub(str_extract(modnames, "&.|&.."),2,-1))
  plot(Q,ICL,xlab="Q1+Q2",ylab="ICL")
}

missforward_explo<-function(models,k1,k2,connectivity,completness,proba_type,estimate_psi,param){
  ##Function to explore all possible sharing combinations of groups from k1 and k2 groups
  cl01<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="&")]]$membership1)
  cl02<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="&")]]$membership2)
  n1=length(unique(cl01))
  n2=length(unique(cl02))
  print(n1)
  print(n2)
  best_one=list()
  if (n1==k1&n2==k2){
    candidates <- mclapply(1:n1, function(j) {##Each group is randomly divided into two groups and the model is run again.
      cl1 <- cl01
      J  <- which(cl1 == j)
      if (length(J) > 1) {
        J1 <- base::sample(J, floor(length(J)/2))
        J2 <- setdiff(J, J1)
        cl1[J1] <- j; 
        cl1[J2] <- n1 + 1
        model=missLBM_VEM(connectivity,completness,n1+1,k2,clustering_indicator(cl1),clustering_indicator(cl02),proba_type,estimate_psi,param)
      }
      else {
        model=models[[paste(as.character(k1+1),as.character(k2),sep="&")]]
      }
      return(model)
    },mc.cores = param$cores)
    best_one[[1]] <- candidates[[which.min(sapply(candidates, function(candidate) candidate$ICL))]]##We take the best of all group splitting in two for the first axis
    
    candidates <- mclapply(1:n2, function(j) {##Same for the second axis
      cl2 <- cl02
      J  <- which(cl2 == j)
      if (length(J) > 1) {
        J1 <- base::sample(J, floor(length(J)/2))
        J2 <- setdiff(J, J1)
        cl2[J1] <- j; 
        cl2[J2] <- n2 + 1
        model=missLBM_VEM(connectivity,completness,k1,k2+1,clustering_indicator(cl01),clustering_indicator(cl2),proba_type,estimate_psi,param)
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



missbackward_explo<-function(models,k1,k2,connectivity,completness,proba_type,estimate_psi,param){
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
    candidates<-mclapply(combn(n1, 2, simplify = FALSE),function(couple){##We explore all the possible combinations of groups for the first axis
      cl_fusion1 <- cl1
      levels(cl_fusion1)[which(levels(cl_fusion1) == paste(couple[1]))] <- paste(couple[2])
      levels(cl_fusion1) <- as.character(1:(n1-1))
      cl_fusion1<-as.numeric(cl_fusion1)
      model=missLBM_VEM(connectivity,completness,n1-1,k2,clustering_indicator(cl_fusion1),clustering_indicator(cl02),proba_type,estimate_psi,param)
      return(model)
    },mc.cores = param$cores)
    best_one[[1]] <- candidates[[which.min(sapply(candidates, function(candidate) candidate$ICL))]]
    
    candidates<-mclapply(combn(n2, 2, simplify = FALSE),function(couple){##We explore all possible combinations of groups for the second axis
      cl_fusion2 <- cl2
      levels(cl_fusion2)[which(levels(cl_fusion2) == paste(couple[1]))] <- paste(couple[2])
      levels(cl_fusion2) <- as.character(1:(n2-1))
      cl_fusion2<-as.numeric(cl_fusion2)
      model=missLBM_VEM(connectivity,completness,k1,n2-1,clustering_indicator(cl01),clustering_indicator(cl_fusion2),proba_type,estimate_psi,param)
      return(model)
    },mc.cores = param$cores)
    best_one[[2]] <- candidates[[which.min(sapply(candidates, function(candidate) candidate$ICL))]]
  }
  else{
    best_one[[1]]=models[[paste(as.character(k1-1),as.character(k2),sep="&")]]##Si le nombre de groupes ne correspond pas à celui exploré, on ne fiat rien.
    best_one[[2]]=models[[paste(as.character(k1),as.character(k2-1),sep="&")]]
  }
  return(best_one)
}



missLBM_update_tau<-function(connectivity,nu,alpha1,alpha2,pi,tau1,tau2){##If the number of groups does not correspond to the number explored, we do not fiat anything.
  barnu=bar_nu(nu,connectivity)
  Q1=length(alpha1)
  Q2=length(alpha2)
  if ((Q1>1)&(Q2>1)){
    newtau1<-connectivity %*% tau2 %*% t(log(pi)) + barnu %*% tau2 %*% t(log(1 - pi)) + nu %*% tau2 %*% t(log(pi))
    newtau1 <- t(apply(sweep(newtau1, 2, log(alpha1), "+"), 1, .softmax))
    newtau2<-t(connectivity) %*% tau1 %*% log(pi) + t(barnu) %*% tau1 %*% log(1 - pi) +t(nu) %*% tau1 %*% log(pi)
    newtau2 <- t(apply(sweep(newtau2, 2, log(alpha2), "+"), 1, .softmax))
  }
  else if ((Q1>1)&(Q2==1)){
    newtau1<-connectivity %*% tau2 %*% t(log(pi)) + barnu %*% tau2 %*% t(log(1 - pi)) + nu %*% tau2 %*% t(log(pi))
    newtau1 <- t(apply(sweep(newtau1, 2, log(alpha1), "+"), 1, .softmax))
    newtau2=tau2
  }
  else if ((Q1==1)&(Q2>1)){
    newtau2<-t(connectivity) %*% tau1 %*% log(pi) + t(barnu) %*% tau1 %*% log(1 - pi) +t(nu) %*% tau1 %*% log(pi)
    newtau2 <- t(apply(sweep(newtau2, 2, log(alpha2), "+"), 1, .softmax))
    newtau1=tau1
  }
  else if ((Q1==1)&(Q2==1)){
    newtau1=tau1
    newtau2=tau2
  }
  return(list(newtau1,newtau2))
}

missLBM_update_nu<-function(samp,connectivity,tau1,tau2,pi){##Function that update nu
  ##In practice we observe that the part tau1%*%.logit(pi)%*%t(tau2) is smalled compared to log(1-samp)
  x=log(1-samp)+tau1%*%.logit(pi)%*%t(tau2)
  nu=.logistic(x)*LBMbar(connectivity)
  return(nu)
}

missLBM_update_pi = function(connectivity,nu,tau1,tau2) { ##We update pi
  N1<-dim(connectivity)[1]
  N2<-dim(connectivity)[2]
  pi    <- check_boundaries((t(tau1)%*%connectivity%*%tau2+t(tau1)%*%nu%*%tau2) / t(tau1)%*%matrix(1,N1,N2)%*%tau2)
  return(pi)
}

missLBM_update_alpha = function(tau) { ##Update of alpha
  alpha <- check_boundaries(colMeans(tau))
  return(alpha)
}

estim_psi<-function(p,psi,Y,nu){
  pbis=(p==1)*-1+p
  res=sum(log(p)*Y)-sum((nu*log(p)*p^psi)/(1-pbis^psi))
  return(res)
}

estim_psi2<-function(p,psi,Y,nu){
  res=sum(Y)/psi-sum((nu*p)/(1-pbis*psi))
  return(res)
}

missLBM_update_psi<-function(connectivity,nu,samp,type_psi){
  if (type_psi==1){
    partial_estim_psi<-function(psi){
      return(estim_psi(samp,psi,connectivity,nu))
    }
    psi=uniroot(partial_estim_psi,interval=c(0.01,4))$root
    samp=samp^psi
  }
  else if (type_psi==2){
    partial_estim_psi2<-function(psi){
      return(estim_psi2(samp,psi,connectivity,nu))
    }
    if (partial_estim_psi2(1)<0){
      psi=uniroot(partial_estim_psi2,interval=c(0.1,1))$root
    }
    else{psi=1}
    samp=psi*samp
  }
  return(list(psi,samp))
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

missLBM_ICL<-function(members1,members2,alpha1,alpha2,pi, connectivity,nu){##Compute ICL from model's data for bernoulli distribution
  Q1=length(alpha1)
  Q2=length(alpha2)
  N1=dim(connectivity)[1]
  N2=dim(connectivity)[2]
  barnu=bar_nu(connectivity,nu)
  logL=sum(members1%*%log(alpha1))+sum(members2%*%log(alpha2))+sum((t(members1)%*%connectivity%*%members2)*log(pi))+sum((t(members1)%*%nu%*%members2)*log(pi))+sum((t(members1)%*%barnu%*%members2)*log(1-pi))
  pena=log(N1)*(Q1-1)/2+log(N1*(N1+1)/2)*Q1*(Q1+1)/4+log(N2)*(Q2-1)/2+log(N2*(N2+1)/2)*Q2*(Q2+1)/4
  return(logL-pena)
}

##Function pour initialiser la probabilité

init_prob<-function(completness,adj,type=2){##Compute the probability of having been sampled for each node from completeness
  ##The completness must be provided in the form of a list of two vectors containing the completeness of axis 1 and axis 2.
  ##or a list containing only the completeness for axis 1 or 2
  ##We have 4 types of probability type1:C1*C2,type2:(C1+C2)/2,type3:C1,type4:C2
  N=dim(adj)
  if (type==1){
    return(completness[[1]]%*%t(completness[[2]]))
  }
  else if (type==2){
    d1=length(completness[[1]])
    d2=length(completness[[2]])
    samp=matrix(0,d1,d2)
    for (i in 1:d1){
      for (j in 1:d2){
        samp[i,j]=(completness[[1]][i]+completness[[2]][j])/2
      }
    }
    return(samp)
  }
  else if (type==3){
    return(completness%x%matrix(1,1,N[2]))
  }
  else if (type==4){
    return(t(completness)%x%matrix(1,N[1],1))
  }
}

missLBM_VEM<-function(connectivity,completness,Q1,Q2,tau1=c(),tau2=c(),proba_type=2,estimate_psi=FALSE,param) {



  ##The completness must be provided in the form of a list of two vectors containing the completeness of axis 1 and 2. 
  ##or a list containing only the completeness for axis 1 or 2
  ##We have 5 types of probability type1:C1*C2,type2:(C1+C2)/2,type3:1-(1-C1)(1-C2),type4:C1,type5:C2
  ##param=list(type=Bernoulli,maxIter=50,fixPointIter=3,threshold=1e-3,trace=TRUE,cores=1)
  ##type in param reprensents the distribution to use in the model
  ##maxIter is the maximum number of iteraction before the algorithlm stop even if the convergence conditions are not reached
  ##fixPointIter is the number of iteraction to use in the fixed point algorithm for tau
  ##threshold is The convergence threshold
  ##trace is the parameter to display some data of the algorithm
  ##cores how many cores to be used
  ##Estimate_psi : False, we do not use this parameter, 1: we use the probability p power psi p^psi, 2: we use the probability psi*p  
  if (param$trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
  ## We initialize the vector of convergences
  delta     <- vector("numeric", param$maxIter)
  ##We initiate the matrix of connection probabilities pi
  pi=matrix(0.5,Q1,Q2)
  ##We calculate the sampling from the completness
  samp=init_prob(completness,connectivity,proba_type)
  sampref=samp
  N=dim(connectivity)
  ##The matrix of the distribution of missing data is initialized. 
  nu=matrix(0,N[1],N[2])
  psi=1
  i <- 0
  cond <- FALSE
  while (!cond) {
    
    i <- i + 1
    
    ## ______________________________________________________
    ## M-step
    #
    # We update the parameters of the LBM alpha and pi
    alpha1=missLBM_update_alpha(tau1)
    alpha2=missLBM_update_alpha(tau2)
    
    new_pi=missLBM_update_pi(connectivity,nu,tau1,tau2)
    
    
    
    ## ______________________________________________________
    ## Variational E-Step
    #
    for (k in seq.int(param$fixPointIter)) {
      tau=missLBM_update_tau(connectivity,nu,alpha1,alpha2,new_pi,tau1,tau2)
      tau1=tau[[1]]
      tau2=tau[[2]]
    }
    ## On met à jour la distribution nu
    nu=missLBM_update_nu(samp,connectivity,tau1,tau2,new_pi)
    
    ## Je mets à jour le paramètre psi
    
    if(estimate_psi!=FALSE){
      est=missLBM_update_psi(connectivity,nu,sampref,estimate_psi)
      psi=est[[1]]
      samp=est[[2]]
    }

    ## We're looking at convergence
    delta[i] <- sqrt(sum((new_pi - pi)^2)) / sqrt(sum((pi)^2))
    cond     <- (i > param$maxIter) |  (delta[i] < param$threshold)
    pi=new_pi
  }
  members1=memberships(tau1)
  members2=memberships(tau2)
  icl=missLBM_ICL(members1,members2,alpha1,alpha2,pi,connectivity,nu)
  ##On rassemble les résultats dans une liste
  res=list()
  res$tau1=tau1
  res$tau2=tau2
  res$pi=pi
  res$alpha1=alpha1
  res$alpha2=alpha2
  res$ICL=icl
  res$membership1=members1##To which group do the species in the first axis belong?
  res$membership2=members2##To which group do the species in the second axis belong?
  res$nu=nu
  res$samp=samp
  res$psi=psi
  return(res)
}

