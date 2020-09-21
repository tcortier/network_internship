require(fields)

miss_LBM_main_VEM<-function(connectivity,completness,proba_type,estimate_psi=FALSE,param,initMethod="hierarchical_clust",maxExplo=1.5,maxGroups=10,reinitialize=FALSE){
  ##Input parameters : 
  ##Adjacency matrix : connectivity
  ##completeness in the form of a list of two vectors or one vector for columns or rows
  ##proba_type : The type of probability to build from completeness : see init_prob in utils_corr_LBM
  ##estimate_psi : whether or not we need to estimate psi
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
      name=unlist(lapply(1:k[2],function(k2){paste(as.character(k[1]),as.character(k2),sep="&")}))##The name of the models in the list
      mods=mclapply(1:k[2],function(k2){##We evaluate all possible combinations with group k1
        print(paste('k={',k[1],',',k2,'}'))
        cl0=clustinit_LBM(connectivity,k[1],k2,initMethod)##We intialize clusters
        tau1=clustering_indicator(cl0[[1]])
        tau2=clustering_indicator(cl0[[2]])
        model<-missLBM_VEM(connectivity,completness,k[1],k2,tau1,tau2,proba_type,estimate_psi,param)
        return(model)
      },mc.cores = param$cores)
      k2max=which.max(sapply(mods, function(mod) mod$ICL))##We see for which number of groups k2 we have the maximum for the ICL
      names(mods)<-name
      models<-c(models,mods)##Estimated models are added to the list of models
      if (models[[paste(as.character(k[1]),as.character(k2max),sep="&")]]$ICL>max){##We update the maximum of the ICL
        whmax=c(k[1],k2max)
        max=models[[paste(as.character(k[1]),as.character(k2max),sep="&")]]$ICL
      }
      LBM_plot(models)##We plot the ICL
    }
    else if (gr==2){##Same thing if we increase the number of k2 groups.
      name=unlist(lapply(1:k[1],function(k1){paste(as.character(k1),as.character(k[2]),sep="&")}))
      mods=mclapply(1:k[1],function(k1){
        print(paste('k={',k1,',',k[2],'}'))
        cl0=clustinit_LBM(connectivity,k1,k[2],initMethod)
        tau1=clustering_indicator(cl0[[1]])
        tau2=clustering_indicator(cl0[[2]])
        model<-missLBM_VEM(connectivity,completness,k1,k[2],tau1,tau2,proba_type,estimate_psi,param)
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
    ##We explore a maximum of 4 groups and then maxExplo times the number of ICL max for k1 and k2
    if ((k[1]<max(4,round((maxExplo*whmax[1])+0.1)))&(k[2]<max(4,round((maxExplo*whmax[2])+0.1)))){##The smallest number of groups is updated if both have not reached the maximum number of participants.
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
  if (reinitialize==TRUE){##If we continue forward and backward exploration
    max2=max
    it<-0
    cond<-TRUE
    print(k)
    while(cond&(it<3)){##We repeat the operation a maximum of 3 times if we have improved the maximum ICL of all our models.
      it<-it+1
      for (k1 in 1:(k[1]-1)){
        for (k2 in 1:(k[2]-1)){
          print(paste('forward','k={',k1,',',k2,'}'))
          model<-missforward_explo(models,k1,k2,connectivity,completness,proba_type,estimate_psi,param)##We explore all combinations and divide a group at random.
          if (model[[1]]$ICL>models[[paste(as.character(k1+1),as.character(k2),sep="&")]]$ICL){
            models[[paste(as.character(k1+1),as.character(k2),sep="&")]]=model[[1]]###If there is a better one we update it in the main list.
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
          model<-missbackward_explo(models,k1,k2,connectivity,completness,proba_type,estimate_psi,param)##All possibilities of grouping from k1 k2 onward are explored.
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
  return(models)##We return the list of models
}
