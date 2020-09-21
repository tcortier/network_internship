#require(Rlab)
#require(LaplacesDemon)
require(CrossClustering)
require(permute)
#require(SpadeR)
library(reshape2)
require(combinat)
require(MixGHD)
#require(ggalluvial)

echantillonnage<-function(adj,P,Z1,Z2,density=list(c(),c()),use_density='unif',sampling_time=c(),sampling_axis=1,freq_interaction='unif',sharpness_of_expo=7){
  ##On réalise un échantillonnage par espèce à l'aide d'une loi de poisson pour chaque interaction possible. 
  ##La probabilité peut également dépendre des densités des espèces. 
  ##La fréquence des interactions peut être uniforme, 'unif', suivre une distribution exponentielle de rareté aléatoirement 'expo',
  ##suivre la probabilité d'interaction des groupes 'groups'...
  n1=dim(adj)[1]
  n2=dim(adj)[2]
  if (use_density=='unif'){
    d1=rep(1,n1)
    d2=rep(1,n2)
    density=list(d1,d2)
  }
  if (freq_interaction=='unif'){
    freq=adj
  }
  if (freq_interaction=='expo'){
    flatadj=c(adj)
    zeros=which(flatadj==1)
    l=length(zeros)
    int_pot=c(1:l)*(sharpness_of_expo/l)
    ind=shuffle(l)
    freq=rep(0,(n1*n2))
    freq[zeros[ind]]=exp(-int_pot)
    freq=matrix(freq,n1,n2)
    freq=freq*adj
  }
  if (freq_interaction=='groups'){
    freq=(Z1%*%P%*%t(Z2))*adj
  }
  if (sampling_axis==1){
    prob=((sampling_time*density[[1]])%*%t(density[[2]]))*freq
  }
  if (sampling_axis==2){
    prob=density[[1]]%*%t(sampling_time*density[[2]])*freq
  }
  sampling=apply(prob,1:2,poiss)
  return(sampling)
}

memb<-function(npc){
  n=sum(npc)
  cl=length(npc)
  Z=matrix(0,n,cl)
  c=c(1,0)
  for(i in 1:cl){
    c[2]=c[2]+npc[i]
    Z[c[1]:c[2],i]=rep(1,npc[i])
    c[1]=c[1]+npc[i]
  }
  return(Z)
}

abundance<-function(N,sharpness=7){
  ind=c(1:N)*(sharpness/N)
  ab=exp(-ind)
  ind=shuffle(N)
  ab=ab[ind]
  return(ab)
}

abundance2<-function(N,sharpness=7){
  ind=c(floor(-N/2+1):floor(N/2))*(sharpness/N)
  ab=exp(-ind)
  ind=shuffle(N)
  ab=ab[ind]
  return(ab)
}

jackniffe2<-function(vect){
  one=length(vect[vect==1])
  two=length(vect[vect==2])
  So=length(vect[vect>0])
  return(max(So,So+2*one-two))
}

compute_comp_jackniffe<-function(freq){
  N=dim(freq)
  estNrow<-c()
  estNcol<-c()
  for (i in 1:N[1]){
    estNrow<-c(estNrow,jackniffe2(freq[i,]))
  }
  for (i in 1:N[2]){
    estNcol<-c(estNcol,jackniffe2(freq[,i]))
  }
  return(list(estNrow,estNcol))
}

poiss=function(x){rpois(1,x)}

ari_comp<-function(mod,ref){
  return(ARI(mod,ref))
}


compute_ari<-function(modelref,model2,iszeros=FALSE,zeros=c(),same_number_of_groups=FALSE){
  ICL1=unlist(lapply(modelref, '[[','ICL'))
  k1=which.max(ICL1)
  k1=names(modelref)[k1]
  if (iszeros==TRUE){
    rcl1=membertoclust(modelref[[k1]]$membership1)[zeros[[1]]]
    ccl1=membertoclust(modelref[[k1]]$membership2)[zeros[[2]]]
    
  }
  else{
    rcl1=membertoclust(modelref[[k1]]$membership1)
    ccl1=membertoclust(modelref[[k1]]$membership2)
  }
  
  
  
  ICL2=unlist(lapply(model2, '[[','ICL'))
  k2=which.max(ICL2)
  if (same_number_of_groups==TRUE){
    rcl2=membertoclust(model2[[k1]]$membership1)
    ccl2=membertoclust(model2[[k1]]$membership2)
  }
  else{
    rcl2=membertoclust(model2[[k2]]$membership1)
    ccl2=membertoclust(model2[[k2]]$membership2)
  }
  arirow=ari_comp(rcl1,rcl2)
  aricol=ari_comp(ccl1,ccl2)
  return(c(arirow,aricol))
}

matrix_distance<-function(modelref,modelcomp){
  ICLref=unlist(lapply(modelref, '[[','ICL'))
  kref=which.max(ICLref)
  kref=names(modelref)[kref]
  matrixref=modelref[[kref]]$pi
  di=dim(matrixref)
  permrow=permn(di[1])
  permcol=permn(di[2])
  matrixcomp=modelcomp[[kref]]$pi
  listdiff<-c()
  for (p1 in permrow){
    for (p2 in permcol){
      diff=sum(sqrt((matrixref-matrixcomp[p1,p2])^2))/length(matrixref)
      listdiff<-c(listdiff,diff)
    }
  }
  return(min(listdiff))
}

