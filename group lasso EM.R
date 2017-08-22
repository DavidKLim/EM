library(caret)
library(mclust)
library(stats)
library(MASS)
library(fpc)
library(permute)
library(flexclust)
library(amap)
library(gplots)
library(glmnet)


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




EM<-function(y,k,alpha=alpha,lambda=lambda,size_factors){

n<-ncol(y)
g<-nrow(y)

vect_y<-as.vector(t(y))
new_y<-rep(vect_y,each=k) # flatten and multiply each count by number of clusters
new_size_factors<-rep(rep(size_factors,each=k),times=g)    ############## corresponding size factors for each entry in trans. data #
gene<-rep(1:g,each=k*n) # gene for each corresponding new_y

#clusts<-matrix(rep(diag(k),times=n*g),byrow=TRUE,ncol=k) # cluster indicators (cell means coding)
design.mat<-matrix(rep(0,times=k^2),ncol=k)
design.mat[,1]<-1
for(c in 1:k){
  design.mat[c,c]<-1
}

clusts<-matrix(rep(t(design.mat),times=n*g),byrow=TRUE,nrow=k*n*g)

# EM
  
  # Initial Clustering
  d<-dist(dist(t(y)))
  #d<-as.dist(1-cor(norm_y, method="spearman"))  ##Spearman correlation distance w/ log transform##
  model<-hclust(d,method="complete")       # hierarchical clustering
  #col<-rep("",times=ncol(y))
  #for(i in 1:length(col)){if(anno$Adeno.Squamous[i]=="adenocarcinoma"){col[i]="red"}else{col[i]="blue"}}
  #heatmap.2(as.matrix(norm_y), Rowv=as.dendrogram(model), Colv=as.dendrogram(model),ColSideColors=col)
  cls<-cutree(model,k=k)
  
  # initialize weights
  wts<-matrix(rep(0,times=k*ncol(y)),nrow=k)
  for(c in 1:k){
    wts[c,]=(cls==c)^2
  }
    
  vect_wts<-rep(as.vector(wts),times=g)
  clust_index<-rep((1:k),times=n*g)
  dat<-cbind(new_y,clusts,clust_index,gene,vect_wts,new_size_factors) # this is k*g*n rows. cols: count, indicator for cl1, cl2, cl3, genes, wts
  
  colnames(dat)[1]<-c("count")
  colnames(dat)[(k+2):ncol(dat)]<-c("clusts","g","weights","size_factors")
  
  
  # initiate finalwts, coefs, maxit, Q, pi
  
  finalwts<-matrix(rep(0,times=k*ncol(y)),nrow=k)
  
  maxit = 100
  
  coefs<-matrix(rep(0,times=g*k),nrow=g)
  pi<-rep(0,times=k)
  Q<-rep(0,times=maxit)
  
  
  
  ########### M / E STEPS #########
  for(a in 1:maxit){
    
    dat[,"weights"]<-rep(as.vector(wts),times=g) # update weights column in dat
    
    # M step
    
    for(j in 1:g){
      dat_g<-dat[dat[,"g"]==j,]                                  # subset just the j'th gene
      
      # USING GLM #
      # for(c in 1:k){
      #   fit<-glm(dat_g[,1]~dat_g[,2:4]-1,weights=dat_g[,ncol(dat_g)])
      #   coefs[j,c]<-log(fit$coefficients[c])
      # }
      
      # GROUP LASSO FROM GLMNET #
      fit<-glmnet(x=dat_g[,2:(k+1)],y=dat_g[,1],family="poisson",alpha=alpha,lambda=lambda,weights=dat_g[,"weights"],offset=log(dat_g[,"size_factors"]),type.multinomial="grouped")   # alpha = 1: lasso
      coefs.glmnet<-coef(fit)[,1]
      coefs[j,]<-coefs.glmnet[1]+coefs.glmnet[2:(k+1)]
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
        l[c,i]<-sum(dpois(y[,i],lambda=exp(coefs[,c])*size_factors[i],log=TRUE))    # posterior log like, include size_factor of subj
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
    
    
  }
  
  
  
  
  
  num_warns=length(warnings())
  
  final_clusters<-rep(0,times=n)
  for(i in 1:n){
    final_clusters[i]<-which.max(finalwts[,i])
  }
  
  m<-rep(0,times=g)
  nondiscriminatory=rep(FALSE,times=g)
  mean_across_clusters<-rowSums(coefs)/ncol(coefs)
  
  for(j in 1:g){
    if(all(abs(coefs[j,]-mean_across_clusters[j])<0.7)){nondiscriminatory[j]=TRUE}     # threshold for nondiscriminatory gene: 1.5 diff from mean across clusters
    for(c in 1:k){
      if(abs(coefs[j,c]-mean_across_clusters[j])>0.7){m[j]=m[j]+1} # nondiscriminatory threshold: away from mean by .7
    }
    if(m[j]==0){m[j]=1}
  }
  
  pred.nondiscriminatory<-mean(nondiscriminatory)
  

  log_L<-sum(apply(log(pi)+l,2,logsumexpc))    # need to check
  
  BIC = -2*log_L + log(n*g)*(sum(m)+(k-1))         # -2log(L) + log(#obs)*(#parameters estimated). minimum = best
  
  result<-list(pi=pi,
               coefs=coefs,
               Q=Q[1:a],
               BIC=BIC,
               nondiscriminatory=nondiscriminatory,
               final_clusters=final_clusters)
  #plot(Q[2:a])   # omit the first one due to instability
  return(result)
  
}
