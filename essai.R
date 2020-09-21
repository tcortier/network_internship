require(blockmodels)
require(igraph)
require(Rlab)
require(LaplacesDemon)
require(CrossClustering)


source("~/Documents/stage_reseaux/codes/correctedlbm/utils_sampling.R")
source("~/Documents/stage_reseaux/codes/correctedlbm/utils_VEM_LBM.R")
source("~/Documents/stage_reseaux/codes/correctedlbm/utils_miss_LBM.R")
source("~/Documents/stage_reseaux/codes/correctedlbm/utils_corr_LBM.R")
source("~/Documents/stage_reseaux/codes/correctedlbm/main_VEM_LBM.R")
source("~/Documents/stage_reseaux/codes/correctedlbm/main_miss_LBM.R")
source("~/Documents/stage_reseaux/codes/correctedlbm/main_corr_LBM.R")

##Generation of an adjacency matrix 

npc1 <- c(10,20,30) # nodes per class
npc2 <- c(20,40,140) # nodes per class
Q <- c(3,3) # classes
n1 <- sum(npc1) # nodes
n2 <- sum(npc2) # nodes
Z1=memb(npc1)
Z2=memb(npc2)
p<-c(0.99,0.8,0.15,0.8,0.2,0.05,0.2,0.1,0.01)
P<-matrix(p,Q[1],Q[2]) 
M<-1*(matrix(runif(n1*n2),n1,n2)<Z1%*%P%*%t(Z2))


abundance_pol=abundance2(n2,sharpness = 6)
abundance_pl=rep(1,n1)
density=list(abundance_pl,abundance_pol)
sampling_time1=runif(n1,min=50,max=150)*0.01
sampl1=echantillonnage(M,P,Z1,Z2,density=density,use_density='defined',sampling_time=sampling_time1,sampling_axis=1,freq_interaction='unif',sharpness_of_expo = 7)
zeroscol=colSums(sampl1)
zeroscol=zeroscol>0
zerosrow=rowSums(sampl1)
zerosrow=zerosrow>0
sampl1=sampl1[,zeroscol]
sampl1=sampl1[zerosrow,]
sampl1bin<-(sampl1>0)*1

compl1<-rowSums(sampl1bin)/rowSums(M[zerosrow,zeroscol])
compl2<-colSums(sampl1bin)/colSums(M[zerosrow,zeroscol])
completude=list(compl1,compl2)


##Parameters for the model
param=list(maxIter=50,fixPointIter=5,threshold=1e-3,trace=TRUE,cores=1,type='bernoulli')

##LBM model

model=LBM_main_VEM(M,param,initMethod="kmeans_clust",maxExplo=1.5,maxGroups=10,reinitialize = FALSE)

##Model miss

modelmiss=miss_LBM_main_VEM(sampl1bin,compl2,proba_type = 4,estimate_psi=FALSE,param=param,initMethod="kmeans_clust",maxExplo=1.5,maxGroups=10,reinitialize=FALSE)

##Model corr

modelcorr=corr_LBM_main_VEM(sampl1bin,compl2,proba_type = 4,estimate_psi=FALSE,param=param,initMethod="kmeans_clust",maxExplo=1.5,maxGroups=10,reinitialize=FALSE)






