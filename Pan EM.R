library(caret)
library(mclust)
library(stats)
library(MASS)
library(fpc)
library(permute)
library(flexclust)
library(amap)
library(gplots)


logsumexpc=function(v){  
  if(any(is.infinite(v))){ 
    stop("infinite value in v\n") 
  } 
  if(length(v)==1){ return(v[1]) }  
  sv = sort(v, decreasing=TRUE) 
  res = sum(exp(sv[-1] - sv[1]))
  lse = sv[1] + log(1+res)
  lse 
}

soft_thresholding=function(alpha,lambda){
  if(abs(alpha)-lambda<0){
    return(0)
  }else{
    return(sign(alpha)*(abs(alpha)-lambda))
  }
}




EM<-function(y, k,
             lambda1=0, lambda2=0, tau=0,
             size_factors=rep(1,times=ncol(y)) ){

n<-ncol(y)
g<-nrow(y)

vect_y<-as.vector(t(y))
new_y<-rep(vect_y,each=k) # flatten and multiply each count by number of clusters
gene<-rep(1:g,each=k*n) # gene for each corresponding new_y
clusts<-matrix(rep(t(diag(k)),times=n*g),byrow=TRUE,ncol=k) # cluster indicators

# EM
  
  # Initial Clustering
  #d<-dist(t(y))                               ##Euclidean distance##
  #d<-dist(t(norm_y))                        ## scaled to account for size factors ##
  d<-as.dist(1-cor(norm_y, method="spearman"))  ##Spearman correlation distance w/ log transform##
  model<-hclust(d,method="complete")       # hierarchical clustering
  #col<-rep("",times=ncol(y))
  #for(i in 1:length(col)){if(anno$Adeno.Squamous[i]=="adenocarcinoma"){col[i]="red"}else{col[i]="blue"}}
  #heatmap.2(as.matrix(norm_y), Rowv=as.dendrogram(model), Colv=as.dendrogram(model),ColSideColors=col)
  cls<-cutree(model,k=k)
  
  #model<-cclust(t(y),k=2) # convex clustering
  #model<-kcca(t(y),k=2)   # K-means
  #model<-kcca(t(y),k=2,family=kccaFamily("kmedians"))   # K-medians
  #model<-kcca(t(y),k=2,family=kccaFamily("angle"))
  #model<-kcca(t(y),k=2,family=kccaFamily("jaccard"))
  #model<-kcca(t(y),k=2,family=kccaFamily("ejaccard"))
  #cls<-predict(model)
  
  #cls<-sample(1:k,n,replace=TRUE) #random initialization
  
  
  # initialize weights
  wts<-matrix(rep(0,times=k*ncol(y)),nrow=k)
  for(c in 1:k){
    wts[c,]=(cls==c)^2
  }
    
  vect_wts<-rep(as.vector(wts),times=g)
  clust_index<-rep((1:k),times=n*g)
  dat<-cbind(new_y,clusts,clust_index,gene,vect_wts) # this is k*g*n rows. cols: count, indicator for cl1, cl2, cl3, genes, wts
  
  colnames(dat)[1]<-c("count")
  colnames(dat)[(k+2):ncol(dat)]<-c("clusts","g","weights")
  
  finalwts<-matrix(rep(0,times=k*ncol(y)),nrow=k)
  
  maxit = 100
  
  coefs<-matrix(rep(0,times=g*k),nrow=g)
  pi<-rep(0,times=k)
  Q<-rep(0,times=maxit)
  theta_list<-list()       # temporary to hold all K x K theta matrices
  
  family=poisson(link="log")       # can specify family here
  #offset=log(size_factors)
  offset=rep(0,times=n)       # no offsets
  
  ########### M / E STEPS #########
  for(a in 1:maxit){
    
    # M step
    
    # estimate parameters by IRLS, gene by gene
    #for(j in 1:g){
    #  for(c in 1:k){
    #  fit = glm(y[j,] ~ 1, weights = wts[c,], family = poisson())   # automatically uses IRLS to find our maximization
    #  coefs[j,c]<-fit$coef
    #  }
    #}
    
    dat[,"weights"]<-rep(as.vector(wts),times=g) # update weights column in dat
    
    # IRWLS:
    
    maxit_IRLS=100
    beta<-rep(0,times=k)
    
    for(j in 1:g){
      if(a==1){
        
        for(c in 1:k){
          beta[c]<-log(mean(as.numeric(y[j,cls==c])))                       # Initialize beta
        }
        
        theta<-matrix(rep(0,times=k^2),nrow=k)
        for(c in 1:k){
          for(cc in 1:k){
            if(theta[c,cc]>=tau){theta[c,cc]<-beta[c]-beta[cc]}             # Initialize theta
            else{theta[c,cc]<-soft_thresholding(beta[c]-beta[cc],lambda2)   # 
            }
          }
        }
        
      }else{
        beta<-coefs[j,]                                                   # Retrieve beta & theta from
        theta<-theta_list[[j]]                                            # previous iteration
      }
      
      temp<-matrix(rep(0,times=maxit_IRLS*k),nrow=maxit_IRLS)    # to test for convergence of IRLS
      dat_g<-dat[dat[,"g"]==j,]                                  # subset just the j'th gene
      
      for(i in 1:maxit_IRLS){
        
        temp[i,]<-beta
        if(a==1 & i==1){
          eta<-matrix(rep(beta,times=n),nrow=n,byrow=TRUE)               # first initialization of eta
        }else if(a>1 & i==1){
          eta<-matrix(rep(beta,times=n),nrow=n,byrow=TRUE) + offset     # Retrieval of eta for IRLS (prev. beta + offset)
        }
        
        for(c in 1:k){
          mu<-exp(eta)
          g_fx<-family$linkinv(eta)
          g_fx_prime<-family$mu.eta(eta)
          var_fx<-family$variance(mu)
          
          dat_gc<-dat_g[dat_g[,"clusts"]==c,]
          
          trans_y<-eta[,c] + dat_gc[,"weights"]*(dat_gc[,"count"]-g_fx[,c])/g_fx_prime[,c] - offset    # subtract size factor from transf. y
        
          
          #beta[c]<-log(glm(dat_gc[,"count"] ~ 1 + offset(log(size_factors,base=2)), weights=dat_gc[,"weights"])$coef)
            
          beta[c]<-((lambda1*((sum(beta)-beta[c])+(sum(theta[c,])-theta[c,c])))+((1/n)*sum(var_fx[,c]*trans_y)))/((lambda1*(k-1))+(1/n)*sum(var_fx[,c]))
          if(beta[c]<(-100)){beta[c]=-100}
          
          eta[,c]<-beta[c] + offset      # add back size factors to eta
        }
        
        # update on theta
        for(c in 1:k){
          for(cc in 1:k){
            if(theta[c,cc]>=tau){theta[c,cc]<-beta[c]-beta[cc]}                      # gTLP from Pan paper
            else{theta[c,cc]<-soft_thresholding(beta[c]-beta[cc],lambda2)}           # 
          }
        }
        
        # break conditions for IRLS
        if(i>1){
          if(sum((temp[i,]-temp[i-1,])^2)<1E-7){
            coefs[j,]<-beta                                 # reached convergence
            theta_list[[j]]<-theta
            break
          }
        }
        if(i==maxit_IRLS){
          coefs[j,]<-beta
          theta_list[[j]]<-theta                           # reached maxit
        }
        
      }
    }

    
    
    
    # update on pi_hat
    for(c in 1:k){
      pi[c]=mean(wts[c,])
      if(pi[c]<1E-6){
        warning(paste("cluster proportion", c, "close to 0"))
        pi[c]=1E-6
      } # lowerbound for pi
      if(pi[c]>(1-1E-6)){
        warning(paste("cluster proportion", c, "close to 1"))
        pi[c]=(1-1E-6)
      } # upperbound for pi
    }
    
    
    # log(f_k(y_i)): summing over all j
    l<-matrix(rep(0,times=k*n),nrow=k)
    for(i in 1:n){
      for(c in 1:k){
        l[c,i]<-sum(dpois(y[,i],lambda=exp(coefs[,c]),log=TRUE))    # posterior log like, include size_factor of subj
      }
    }
    
    # store and check stopping criterion
    pt1<-(log(pi)%*%rowSums(wts))
    pt2<-sum(wts*l)
    Q[a]<-pt1+pt2
    
  
    
    
    
    if(a>10){if(abs(Q[a]-Q[a-10])<1E-5) {
      finalwts<-wts
      break
    }}
    
    # E step
    
    # update on weights
    logdenom = apply(log(pi) + l, 2,logsumexpc)
    for(c in 1:k){
      wts[c,]<-exp(log(pi[c])+l[c,]-logdenom)
    }
    
    if(a==maxit){finalwts<-wts}
    # print(pi) # print estimated cluster proportions
    
    # if(any(rowSums(wts)==0)){
    #   print(paste("Empty cluster when K =",k,". Choose smaller K"))
    #   break
    # }
    
  }
  
  
  
  
  
num_warns=length(warnings())
  
  final_clusters<-rep(0,times=n)
  for(i in 1:n){
    final_clusters[i]<-which.max(finalwts[,i])
  }
  
  m<-rep(0,times=g)
  nondiscriminatory=rep(FALSE,times=g)
  #mean_across_clusters<-rowSums(coefs)/ncol(coefs)
  
  for(j in 1:g){
    #if(all(abs(coefs[j,]-mean_across_clusters[j])<0.7)){nondiscriminatory[j]=TRUE}     # threshold for nondiscriminatory gene: 1.5 diff from mean across clusters
    #for(c in 1:k){
    #  if(abs(coefs[j,c]-mean_across_clusters[j])>0.7){m[j]=m[j]+1} # nondiscriminatory threshold: away from mean by 7
    #}
    m[j] <- sum(theta_list[[j]][1,]!=0) + 1         # of parameters estimated
    if(m[j]==1){nondiscriminatory[j]=TRUE}
  }
  
  pred.nondiscriminatory<-mean(nondiscriminatory)
  

  log_L<-sum(apply(log(pi) + l, 2, logsumexpc))
  
  BIC=-2*log_L+log(n*g)*(sum(m)+(k-1))         # -2log(L) + log(#obs)*(#parameters estimated). minimum = best. g*k: total params, sum(m): total # of discriminatory genes
  
  result<-list(pi=pi,
               coefs=coefs,
               Q=Q[1:a],
               BIC=BIC,
               nondiscriminatory=nondiscriminatory,
               final_clusters=final_clusters)
  #plot(Q[2:a])   # omit the first one due to instability
  return(result)
  
}
