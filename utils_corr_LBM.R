require(fields)
require(parallel)
require(Rcpp)

##Here set the directory of the C++ code 
setwd("~/Documents/stage_reseaux/codes/correctedlbm/")
sourceCpp("corrLBM.cpp")


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

check_boundaries_2 <- function(x, zero = .Machine$double.eps) {
  x[is.nan(x)] <- zero
  x[x <     zero] <-     zero
  x
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

corrforward_explo<-function(models,k1,k2,connectivity,completness,proba_type,estimate_psi,param){
  ##Function to explore all possible sharing combinations of groups from k1 and k2 groups
  cl01<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="&")]]$membership1)
  cl02<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="&")]]$membership2)
  n1=length(unique(cl01))
  n2=length(unique(cl02))
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
        model=corrLBM_VEM(connectivity,completness,n1+1,k2,clustering_indicator(cl1),clustering_indicator(cl02),proba_type,estimate_psi,param)
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
        model=corrLBM_VEM(connectivity,completness,k1,k2+1,clustering_indicator(cl01),clustering_indicator(cl2),proba_type,estimate_psi,param)
      }
      else {
        model=models[[paste(as.character(k1),as.character(k2+1),sep="&")]]
      }
      return(model)
    },mc.cores = param$cores)
    best_one[[2]] <- candidates[[which.min(sapply(candidates, function(candidate) candidate$ICL))]]
  }
  else{
    best_one[[1]]=models[[paste(as.character(k1+1),as.character(k2),sep="&")]]##If you don't have the right number of groups you don't do anything.
    best_one[[2]]=models[[paste(as.character(k1),as.character(k2+1),sep="&")]]
  }
  return(best_one)
}



corrbackward_explo<-function(models,k1,k2,connectivity,completness,proba_type,estimate_psi,param){
  cl01<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="&")]]$membership1)
  cl02<-membertoclust(models[[paste(as.character(k1),as.character(k2),sep="&")]]$membership2)
  cl1 <- factor(cl01)
  cl2 <- factor(cl02)
  n1=nlevels(cl1)
  n2=nlevels(cl2)
  best_one=list()
  if (n1==k1&n2==k2){
    candidates<-mclapply(combn(n1, 2, simplify = FALSE),function(couple){##We explore all the possible combinations of groups for the first axis
      cl_fusion1 <- cl1
      levels(cl_fusion1)[which(levels(cl_fusion1) == paste(couple[1]))] <- paste(couple[2])
      levels(cl_fusion1) <- as.character(1:(n1-1))
      cl_fusion1<-as.numeric(cl_fusion1)
      model=corrLBM_VEM(connectivity,completness,n1-1,k2,clustering_indicator(cl_fusion1),clustering_indicator(cl02),proba_type,estimate_psi,param)
      return(model)
    },mc.cores = param$cores)
    best_one[[1]] <- candidates[[which.min(sapply(candidates, function(candidate) candidate$ICL))]]
    
    candidates<-mclapply(combn(n2, 2, simplify = FALSE),function(couple){##We explore all possible combinations of groups for the second axis
      cl_fusion2 <- cl2
      levels(cl_fusion2)[which(levels(cl_fusion2) == paste(couple[1]))] <- paste(couple[2])
      levels(cl_fusion2) <- as.character(1:(n2-1))
      cl_fusion2<-as.numeric(cl_fusion2)
      model=corrLBM_VEM(connectivity,completness,k1,n2-1,clustering_indicator(cl01),clustering_indicator(cl_fusion2),proba_type,estimate_psi,param)
      return(model)
    },mc.cores = param$cores)
    best_one[[2]] <- candidates[[which.min(sapply(candidates, function(candidate) candidate$ICL))]]
  }
  else{
    best_one[[1]]=models[[paste(as.character(k1-1),as.character(k2),sep="&")]]##If the number of groups does not correspond to the number explored, we do not fiat anything.
    best_one[[2]]=models[[paste(as.character(k1),as.character(k2-1),sep="&")]]
  }
  return(best_one)
}



corrLBM_update_tau_poisson<-function(connectivity,sampl,alpha1,alpha2,lambda,tau1,tau2){##Update of tau1 and tau2 taking into account sampling for a poisson distribution
  N1=dim(barconnectivity)[1]
  N2=dim(barconnectivity)[2]
  Q1=length(alpha1)
  Q2=length(alpha2)
  
  if ((Q1>1)&(Q2>1)){
    newtau1=matrix(0,N1,Q1)
    newtau2=matrix(0,N2,Q2)
    for (i in 1:N1){
      for (q in 1:Q1){
        s=log(alpha1[q])
        for (j in 1:N2){
          for (l in 1:Q2){
            s=s+tau2[j,l]*connectivity[i,j]*log(sampl[i,j]*lambda[q,l])-tau2[j,l]*sampl[i,j]*lambda[q,l]-tau2[j,l]*log(factorial(connectivity[i,j]))
          }
        }
        newtau1[i,q]=s
      }
    }
    for (j in 1:N2){
      for (l in 1:Q2){
        s=log(alpha2[l])
        for (i in 1:N1){
          for (q in 1:Q1){
            s=s+tau1[i,q]*connectivity[i,j]*log(sampl[i,j]*lambda[q,l])-tau1[i,q]*sampl[i,j]*lambda[q,l]-tau1[i,q]*log(factorial(connectivity[i,j]))
          }
        }
        newtau2[j,l]=s
      }
    }
    
    newtau1<-t(apply(newtau1, 1, .softmax))
    newtau2<-t(apply(newtau2, 1, .softmax))
  }
  else if ((Q1>1)&(Q2==1)){
    newtau1=matrix(0,N1,Q1)
    for (i in 1:N1){
      for (q in 1:Q1){
        s=log(alpha1[q])
        for (j in 1:N2){
          for (l in 1:Q2){
            s=s+tau2[j,l]*connectivity[i,j]*log(sampl[i,j]*lambda[q,l])-tau2[j,l]*sampl[i,j]*lambda[q,l]-tau2[j,l]*log(factorial(connectivity[i,j]))
          }
        }
        newtau1[i,q]=s
      }
    }
    newtau1<-t(apply(newtau1, 1, .softmax))
    newtau2=tau2
  }
  else if ((Q1==1)&(Q2>1)){
    newtau2=matrix(0,N2,Q2)
    for (j in 1:N2){
      for (l in 1:Q2){
        s=log(alpha2[l])
        for (i in 1:N1){
          for (q in 1:Q1){
            s=s+tau1[i,q]*connectivity[i,j]*log(sampl[i,j]*lambda[q,l])-tau1[i,q]*sampl[i,j]*lambda[q,l]-tau1[i,q]*log(factorial(connectivity[i,j]))
          }
        }
        newtau2[j,l]=s
      }
    }
    newtau2<-t(apply(newtau2, 1, .softmax))
    newtau1=tau1
  }
  else if ((Q1==1)&(Q2==1)){
    newtau1=tau1
    newtau2=tau2
  }
  return(list(newtau1,newtau2))
}


corrLBM_update_tau<-function(connectivity,sampl,alpha1,alpha2,pi,tau1,tau2){##Update of tau1 and tau2 taking into account sampling for a binomiale distribution
  barconnectivity=LBMbar(connectivity)
  N1=dim(barconnectivity)[1]
  N2=dim(barconnectivity)[2]
  Q1=length(alpha1)
  Q2=length(alpha2)
  
  if ((Q1>1)&(Q2>1)){
    newtau1=matrix(0,N1,Q1)
    newtau2=matrix(0,N2,Q2)
    for (i in 1:N1){
      for (q in 1:Q1){
        s=log(alpha1[q])
        for (j in 1:N2){
          for (l in 1:Q2){
            s=s+tau2[j,l]*connectivity[i,j]*log(sampl[i,j]*pi[q,l])+tau2[j,l]*barconnectivity[i,j]*log(1-sampl[i,j]*pi[q,l])
          }
        }
        newtau1[i,q]=s
      }
    }
    for (j in 1:N2){
      for (l in 1:Q2){
        s=log(alpha2[l])
        for (i in 1:N1){
          for (q in 1:Q1){
            s=s+tau1[i,q]*connectivity[i,j]*log(sampl[i,j]*pi[q,l])+tau1[i,q]*barconnectivity[i,j]*log(1-sampl[i,j]*pi[q,l])
          }
        }
        newtau2[j,l]=s
      }
    }
    
    newtau1<-t(apply(newtau1, 1, .softmax))
    newtau2<-t(apply(newtau2, 1, .softmax))
  }
  else if ((Q1>1)&(Q2==1)){
    newtau1=matrix(0,N1,Q1)
    for (i in 1:N1){
      for (q in 1:Q1){
        s=log(alpha1[q])
        for (j in 1:N2){
          for (l in 1:Q2){
            s=s+tau2[j,l]*connectivity[i,j]*log(sampl[i,j]*pi[q,l])+tau2[j,l]*barconnectivity[i,j]*log(1-sampl[i,j]*pi[q,l])
          }
        }
        newtau1[i,q]=s
      }
    }
    newtau1<-t(apply(newtau1, 1, .softmax))
    newtau2=tau2
  }
  else if ((Q1==1)&(Q2>1)){
    newtau2=matrix(0,N2,Q2)
    for (j in 1:N2){
      for (l in 1:Q2){
        s=log(alpha2[l])
        for (i in 1:N1){
          for (q in 1:Q1){
            s=s+tau1[i,q]*connectivity[i,j]*log(sampl[i,j]*pi[q,l])+tau1[i,q]*barconnectivity[i,j]*log(1-sampl[i,j]*pi[q,l])
          }
        }
        newtau2[j,l]=s
      }
    }
    newtau2<-t(apply(newtau2, 1, .softmax))
    newtau1=tau1
  }
  else if ((Q1==1)&(Q2==1)){
    newtau1=tau1
    newtau2=tau2
  }
  return(list(newtau1,newtau2))
}


corrLBM_update_pi = function(connectivity,sampl,tau1,tau2) { ##We update pi for a binomial distribution
  N1<-dim(connectivity)[1]
  N2<-dim(connectivity)[2]
  pi    <- check_boundaries(t(tau1)%*%connectivity%*%tau2 / t(tau1)%*%sampl%*%tau2)
  return(pi)
}

corrLBM_update_lambda = function(connectivity,sampl,tau1,tau2) { ##Update lambda for a poisson distribution
  N1<-dim(connectivity)[1]
  N2<-dim(connectivity)[2]
  lambda    <- check_boundaries_2(t(tau1)%*%connectivity%*%tau2 / t(tau1)%*%sampl%*%tau2)
  return(lambda)
}

estimate_psi<-function(tau1,tau2,pi,connectivity,sampl,psi){
  res<-0
  Q1=dim(tau1)[2]
  Q2=dim(tau2)[2]
  N1=dim(connectivity)[1]
  N2=dim(connectivity)[2]
  for (i in 1:N1){
    for (q in 1:Q1){
      for (j in 1:N2){
        for (l in 1:Q2){
          res=res+tau1[i,q]*tau2[j,l]*((pi[q,l]^connectivity[i,j])*connectivity[i,j]*log(sampl[i,j])*(sampl[i,j]^(psi*connectivity[i,j]))-(1-connectivity[i,j])*log(sampl[i,j])*(sampl[i,j]^psi)*pi[q,l]*(1-(sampl[i,j]^psi)*pi[q,l])^(-connectivity[i,j]))  
        }
      }
    }
  }
  return(res)
}

corrLBM_update_psi = function(connectivity,sampl,tau1,tau2,pi) { ##Update of psi
  partial_estimate_psi<-function(psi){
    return(estimate_psi(tau1,tau2,pi,connectivity,sampl,psi))
  }
  if (partial_estimate_psi(4)<0){
    psi=uniroot(partial_estimate_psi,interval=c(0.1,4))$root
  }
  else{psi=1}
  return(list(psi,sampl^psi))
}

corrLBM_update_alpha = function(tau) { ##Update alpha
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

corrLBM_ICL<-function(members1,members2,alpha1,alpha2,pi, connectivity,sampl){##Compute ICL from model's data for bernoulli distribution
  Q1=length(alpha1)
  Q2=length(alpha2)
  N1=dim(connectivity)[1]
  N2=dim(connectivity)[2]
  logL=corr_compute_logL(members1,members2,pi, connectivity,sampl)
  logL=logL+sum(members1%*%log(alpha1))+sum(members2%*%log(alpha2))
  pena=log(N1)*(Q1-1)/2+log(N1*(N1+1)/2)*Q1*(Q1+1)/4+log(N2)*(Q2-1)/2+log(N2*(N2+1)/2)*Q2*(Q2+1)/4
  
  return(logL-pena)
}


corrLBM_ICL_poisson<-function(members1,members2,alpha1,alpha2,lambda, connectivity,sampl){##Compute ICL from model's data for poisson distribution
  Q1=length(alpha1)
  Q2=length(alpha2)
  N1=dim(connectivity)[1]
  N2=dim(connectivity)[2]
  logL=corr_compute_logL_poisson(members1,members2,lambda, connectivity,sampl)
  logL=logL+sum(members1%*%log(alpha1))+sum(members2%*%log(alpha2))
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


corrLBM_VEM<-function(connectivity,completness,Q1,Q2,tau1=c(),tau2=c(),proba_type=2,estimate_psi=FALSE,param) {
  if (param$type=='bernoulli'){
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
    
    if (param$trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
    ## We initialize the vector of convergences
    delta     <- vector("numeric", param$maxIter)
    ##We initiate the matrix of connection probabilities pi
    pi=matrix(0.5,Q1,Q2)
    psi=1
    ##We calculate the sampling from the completness
    sampl=init_prob(completness,connectivity,proba_type)
    sampl_ref=sampl
    N=dim(connectivity)
    i <- 0
    cond <- FALSE
    while (!cond) {
      
      i <- i + 1
      
      ## ______________________________________________________
      ## M-step
      #
      # We update the parameters of the LBM alpha and pi
      alpha1=corrLBM_update_alpha(tau1)
      alpha2=corrLBM_update_alpha(tau2)
      
      new_pi=corrLBM_update_pi(connectivity,sampl,tau1,tau2)
      
      
      
      ## ______________________________________________________
      ## Variational E-Step
      #
      
      for (k in seq.int(param$fixPointIter)) {
        if (dim(tau1)[2]>1){
          new_tau1=corr_update_tau1(connectivity,sampl,alpha1,new_pi,tau2)
        }
        else{
          new_tau1=tau1
        }
        if (dim(tau2)[2]>1){
          tau2=corr_update_tau2(connectivity,sampl,alpha2,new_pi,tau1)
        }
        tau1=new_tau1
      }
      
      if (estimate_psi==TRUE){
        res=corrLBM_update_psi(connectivity,sampl_ref,tau1,tau2,new_pi)
        psi=res[[1]]
        sampl=res[[2]]
      }
      
      
      ## We're looking at convergence
      delta[i] <- sqrt(sum((new_pi - pi)^2)) / sqrt(sum((pi)^2))
      cond     <- (i > param$maxIter) |  (delta[i] < param$threshold)
      pi=new_pi
    }
    members1=memberships(tau1)
    members2=memberships(tau2)
    icl=corrLBM_ICL(members1,members2,alpha1,alpha2,pi,connectivity,sampl)
    ##We gather the results in a list
    res=list()
    res$tau1=tau1
    res$tau2=tau2
    res$pi=pi
    res$alpha1=alpha1
    res$alpha2=alpha2
    res$ICL=icl
    res$membership1=members1##To which group do the species in the first axis belong?
    res$membership2=members2##To which group do the species in the second axis belong?
    res$samp=sampl
    res$psi=psi
    return(res)
  }
  else if(param$type=='poisson'){
    if (param$trace) cat("\n Adjusting Variational EM for Stochastic Block Model\n")
    delta     <- vector("numeric", param$maxIter)
    lambda=matrix(0,Q1,Q2)
    sampl=init_prob(completness,connectivity,proba_type)
    sampl_ref=sampl
    N=dim(connectivity)
    psi=1

    i <- 0
    cond <- FALSE
    while (!cond) {
      
      i <- i + 1
      
      ## ______________________________________________________
      ## M-step
      #
      # We update the parameters of the LBM alpha and lambda
      alpha1=corrLBM_update_alpha(tau1)
      alpha2=corrLBM_update_alpha(tau2)
      
      new_lambda=corrLBM_update_lambda(connectivity,sampl,tau1,tau2)
      
      
      
      ## ______________________________________________________
      ## Variational E-Step
      #
      for (k in seq.int(param$fixPointIter)) {
        if (dim(tau1)[2]>1){
          new_tau1=corr_update_tau1_poisson(connectivity,sampl,alpha1,new_lambda,tau2)
        }
        else{
          new_tau1=tau1
        }
        if (dim(tau2)[2]>1){
          tau2=corr_update_tau2_poisson(connectivity,sampl,alpha2,new_lambda,tau1)
        }
        tau1=new_tau1
      }
      
      if (estimate_psi==TRUE){
        res=corrLBM_update_psi(connectivity,sampl_ref,tau1,tau2,new_lambda)
        psi=res[[1]]
        sampl=res[[2]]
      }
      
      delta[i] <- sqrt(sum((new_lambda - lambda)^2)) / sqrt(sum((lambda)^2))
      cond     <- (i > param$maxIter) |  (delta[i] < param$threshold)
      lambda=new_lambda
    }
    members1=memberships(tau1)
    members2=memberships(tau2)
    icl=corrLBM_ICL_poisson(members1,members2,alpha1,alpha2,lambda,connectivity,sampl)

    res=list()
    res$tau1=tau1
    res$tau2=tau2
    res$lambda=lambda
    res$alpha1=alpha1
    res$alpha2=alpha2
    res$ICL=icl
    res$membership1=members1
    res$membership2=members2
    res$samp=sampl
    res$psi=psi
    return(res)
  }
}


