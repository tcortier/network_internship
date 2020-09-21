
require(fields)
require(parallel)

LBM_main_VEM<-function(connectivity,param,initMethod="hierarchical_clust",maxExplo=1.5,maxGroups=10,reinitialize=TRUE){
  ##Input parameters : 
  ##Adjacency matrix : connectivity
  ##Param the parameters of the VEM algorithm
  ##initMethod : the method of group initalization "hierarchical_clust", "spectral_clust","kmeans_clust"
  ##MaxExplo : How many groups we explore after the maximum of the ICL
  ##maxGroups : The maximum number of groups to be explored
  ##reinitialize :  whether to re-explore by merging the number of groups or by splitting groups in half (False or True)
  #------
  k<-c(1,1)##Number of groups for columns and rows
  models<-list()##The list containing all the information about the models
  max=-1e20##The maximum of ICL
  whmax=c(0,0)##The number of groups for which we have the ICL max.
  gr=1##If we increase k for rows or columns (1 or 2)
  cond=TRUE##The condition as long as one continues to explore combinations of groups
  while(cond){
    if (gr==1){
      name=unlist(lapply(1:k[2],function(k2){paste(as.character(k[1]),as.character(k2),sep="&")}))
      mods=mclapply(1:k[2],function(k2){
        print(paste('k={',k[1],',',k2,'}'))
        cl0=clustinit_LBM(connectivity,k[1],k2,initMethod)
        tau1=clustering_indicator(cl0[[1]])
        tau2=clustering_indicator(cl0[[2]])
        model<-LBM_VEM(connectivity,k[1],k2,tau1,tau2,param)
        return(model)
      },mc.cores = param$cores)
      k2max=which.max(sapply(mods, function(mod) mod$ICL))
      names(mods)<-name
      models<-c(models,mods)
      if (models[[paste(as.character(k[1]),as.character(k2max),sep="&")]]$ICL>max){
        whmax=c(k[1],k2max)
        max=models[[paste(as.character(k[1]),as.character(k2max),sep="&")]]$ICL
      }
      LBM_plot(models)
    }
    else if (gr==2){
      name=unlist(lapply(1:k[1],function(k1){paste(as.character(k1),as.character(k[2]),sep="&")}))
      mods=mclapply(1:k[1],function(k1){
        print(paste('k={',k1,',',k[2],'}'))
        cl0=clustinit_LBM(connectivity,k1,k[2],initMethod)
        tau1=clustering_indicator(cl0[[1]])
        tau2=clustering_indicator(cl0[[2]])
        model<-LBM_VEM(connectivity,k1,k[2],tau1,tau2,param)
        return(model)
      },mc.cores = param$cores)
      k1max=which.max(sapply(mods, function(mod) mod$ICL))
      names(mods)<-name
      models<-c(models,mods)
      if (models[[paste(as.character(k1max),as.character(k[2]),sep="&")]]$ICL>max){
        whmax=c(k1max,k[2])
        max=models[[paste(as.character(k1max),as.character(k[2]),sep="&")]]$ICL
      }
      LBM_plot(models)
    }
    cond=((k[1]<4|k[1]<round((maxExplo*whmax[1])+0.1)|k[2]<4|k[2]<round((maxExplo*whmax[2])+0.1))&(k[1]<maxGroups)&(k[2]<maxGroups))
    if ((k[1]<max(4,round((maxExplo*whmax[1])+0.1)))&(k[2]<max(4,round((maxExplo*whmax[2])+0.1)))){
      gr=which.min(k)
      k[which.min(k)]<-k[which.min(k)]+1
    }
    else if((k[1]>=max(4,round((maxExplo*whmax[1])+0.1)))&(k[2]<max(4,round((maxExplo*whmax[2])+0.1)))){
      k[2]<-k[2]+1
      gr=2
    }
    else if((k[1]<max(4,round((maxExplo*whmax[1])+0.1)))&(k[2]>=max(4,round((maxExplo*whmax[2])+0.1)))){
      k[1]<-k[1]+1
      gr=1
    }
  }
  if (reinitialize==TRUE){
    max2=max
    it<-0
    cond<-TRUE
    print(k)
    while(cond&(it<3)){
      it<-it+1
      for (k1 in 1:(k[1]-1)){
        for (k2 in 1:(k[2]-1)){
          print(paste('forward','k={',k1,',',k2,'}'))
          model<-forward_explo(models,k1,k2,connectivity,param)
          if (model[[1]]$ICL>models[[paste(as.character(k1+1),as.character(k2),sep="&")]]$ICL){
            models[[paste(as.character(k1+1),as.character(k2),sep="&")]]=model[[1]]
            if (models[[paste(as.character(k1+1),as.character(k2),sep="&")]]$ICL>max2){
              max2=models[[paste(as.character(k1+1),as.character(k2),sep="&")]]$ICL
            }
          }
          
          if (model[[2]]$ICL>models[[paste(as.character(k1),as.character(k2+1),sep="&")]]$ICL){
            models[[paste(as.character(k1),as.character(k2+1),sep="&")]]=model[[2]]
            if (models[[paste(as.character(k1),as.character(k2+1),sep="&")]]$ICL>max2){
              max2=models[[paste(as.character(k1),as.character(k2+1),sep="&")]]$ICL
            }
          }
          LBM_plot(models)
        }
      }
      for (k1 in c(k[1]:3)){
        for (k2 in c(k[2]:3)){
          print(paste('backward','k={',k1,',',k2,'}'))
          model<-backward_explo(models,k1,k2,connectivity,param)
          if (model[[1]]$ICL>models[[paste(as.character(k1-1),as.character(k2),sep="&")]]$ICL){
            models[[paste(as.character(k1-1),as.character(k2),sep="&")]]=model[[1]]
            if (models[[paste(as.character(k1-1),as.character(k2),sep="&")]]$ICL>max2){
              max2=models[[paste(as.character(k1-1),as.character(k2),sep="&")]]$ICL
            }
          }
          if (model[[2]]$ICL>models[[paste(as.character(k1),as.character(k2-1),sep="&")]]$ICL){
            models[[paste(as.character(k1),as.character(k2-1),sep="&")]]=model[[2]]
            if (models[[paste(as.character(k1),as.character(k2-1),sep="&")]]$ICL>max2){
              max2=models[[paste(as.character(k1),as.character(k2-1),sep="&")]]$ICL
            }
          }
          LBM_plot(models)
        }
      }
      if (max2>max){
        max=max2
      }
      else{
        cond=FALSE
      }
    }
  }
  return(models)
}


