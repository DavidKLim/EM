library(caret)
library(mclust)
library(stats)
library(MASS)
library(fpc)
library(permute)


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


EM<-function(n,k,g,init_pi,b,alpha,lambda1,lambda2){

y<-matrix(rep(0,times=g*n),nrow=g)      # initialize count matrix gxn #

# Prepare new flattened data

z = rmultinom(n,1,init_pi)
while(any(rowSums(z)==0)){z=rmultinom(n,1,init_pi)}   # makes sure that no one cluster simulated @ 0 membership

for(j in 1:g){
  for(c in 1:k){
    y[j,z[c,]==1] = rpois(sum(z[c,]==1), lambda = exp(b[j,c]))
  }
}
true_clusters<-rep(0,times=n)
for(i in 1:n){
  true_clusters[i]<-which(z[,i]==1)
}

print(true_clusters)     # print the cluster memberships being simulated

vect_y<-as.vector(t(y))
new_y<-rep(vect_y,each=k) # flatten and multiply each count by number of clusters
gene<-rep(1:g,each=k*n) # gene for each corresponding new_y
clusts<-matrix(rep(diag(k),times=n*g),byrow=TRUE,ncol=k) # cluster indicators

# EM
  
  # Clustering
  d<-dist(t(y))
  model<-hclust(d)
  cls<-cutree(model,k=k)
  #cls<-sample(1:k,n,replace=TRUE) #random initialization
  
  # choosing the permutation of cluster indices that best fit the true cluster indices (for summarizing simulations)
  all_perms=allPerms(1:k)
  all_clusts=list()
  temp_clust<-rep(0,times=n)
  
  for(ii in 1:nrow(all_perms)){
    for(i in 1:n){
      temp_clust[i]<-all_perms[ii,cls[i]]
    }
    all_clusts[[ii]]<-temp_clust
  }
  
  all_clusts[[nrow(all_perms)+1]]<-cls     # contains all permutations of cluster indices
  match_index<-rep(0,times=nrow(all_perms)+1)
  
  for(ii in 1:(nrow(all_perms)+1)){
    match_index[ii]<-mean(true_clusters==all_clusts[[ii]])     # compares each permutation to true --> % of match
  }

  cls<-all_clusts[[which.max(match_index)]]
  
  
  # initialize weights
  wts<-matrix(rep(0,times=k*ncol(y)),nrow=k)
  for(c in 1:k){
    wts[c,]=(cls==c)^2
  }
  
  vect_wts<-rep(as.vector(wts),times=g)
  clust_index<-rep(c(1,2,3),times=n*g)
  dat<-cbind(new_y,clusts,clust_index,gene,vect_wts) # this is k*g*n rows. cols: count, indicator for cl1, cl2, cl3, genes, wts
  
  colnames(dat)<-c("count","cl1","cl2","cl3","clusts","g","weights") # breaks down k>3
  
  finalwts<-matrix(rep(0,times=k*ncol(y)),nrow=k)
  
  maxit = 100
  
  coefs<-matrix(rep(0,times=g*k),nrow=g)
  pi<-rep(0,times=k)
  Q<-rep(0,times=maxit)
  
  
  
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
    
    maxit_IRLS=1000
    for(j in 1:g){
      if(a>1) {eta<-coefs[j,]} else {eta<-rep(0.5,times=k)}
      temp<-matrix(rep(0,times=maxit_IRLS*k),nrow=maxit_IRLS)
      dat_g<-dat[dat[,"g"]==j,]      # subset just the j'th gene
      
      # USING GLM #
      # for(c in 1:k){
      #   fit<-glm(dat_g[,1]~dat_g[,2:4]-1,weights=dat_g[,ncol(dat_g)])
      #   coefs[j,c]<-log(fit$coefficients[c])
      # }
      
      theta<-matrix(rep(0,times=k^2),nrow=k) #initialize theta (beta_k - beta_l)
      
      for(i in 1:maxit_IRLS){
        
        mu=exp(eta)
        temp[i,]<-eta
        
        for(c in 1:k){
          dat_gc<-dat_g[dat_g[,"clusts"]==c,]
          X<-dat_gc[,c+1]
          W<-diag(rep(mu[c],times=n))
          trans_y<-rep(0,times=n)
          trans_y<-eta[c]+dat_gc[,"weights"]*(dat_gc[,"count"]-mu[c])/mu[c]
        
          #eta_update<-ginv(t(X) %*% W %*% X) %*% t(X) %*% W %*% trans_y    #no penalty
          #eta_update<-(lambda*(sum(eta)-eta[c])-trans_y %*% W %*% X)/(2*lambda-n*mu[c]) # penalty all wrong
          
          eta_update1<-lambda1*(sum(eta)-eta[c]+sum(theta[c,]-theta[c,c]))      #second term in numerator
          eta_update2<-(1/n)*sum(mu[c]*trans_y)                    #first term in numerator
          eta_update3<-lambda1*(k-1)                                            #first term in denominator
          eta_update4<-mu[c]*nrow(dat_gc)/n                                           #second term in denominator
            
          eta_update<-(eta_update1+eta_update2)/(eta_update3+eta_update4)
          
          eta[c]<-eta_update
        }
        
        # update on theta
        for(c in 1:k){
          for(cc in 1:k){
            #theta[c,cc]<-soft_thresholding(eta[c]-eta[cc],lambda2/lambda1)   # original Wei Pan Lasso
            theta[c,cc]<-soft_thresholding(eta[c]-eta[cc],alpha*lambda2)/(1+(lambda2/lambda1)*(1-alpha))   # elastic net
          }
        }
        
        # break condition for IRLS
        if(i>1){
          if(sum((temp[i,]-temp[i-1,])^2)<1E-7){
            coefs[j,]<-eta
            break
          }
        }
        if(i==maxit_IRLS){coefs[j,]<-eta}
      }
    }
        ##########
        
        
        #for (c in 1:k){
        #  nz<-eta[c]+wts[c,]*(y[j,]-mu[c])/mu[c]  # weighted observations (?) based on current weights w_jk
          # w[c,]<-rep(mu[c],times=length(nz)) # iterative weights for IRLS (not used)
          #m=0
          #for(d in 1:k){
          #  if(d != c){if(eta[c] != eta[d]) {m=m+abs(eta[c]-eta[d])}}
          #}
        #  m=0
        #  for(f in 1:k){
        #    if(abs(eta[c]-eta[f])>0.05) {m=m+1} # threshold of cluster means equal for penalty: 0.05
        #  }
        #  eta[c]<-mean(nz)       # gives negative value for last one (c=3)
        #  mu[c]<-exp(eta[c])
        #  temp[i,c]<-eta[c]
        #}
        #if(i>1){
        #  if(sum((temp[i,]-temp[i-1,])^2)<1E-7){
        #    coefs[j,]<-eta
        #    break
        #  }
        #}
        #if(i==maxit_IRLS){coefs[j,]<-eta}
      #}
    #}
    
    
    
    # update on pi_hat
    for(c in 1:k)
    {
      pi[c]=mean(wts[c,])
    }
    
    
    # log(f_k(y_i))
    l<-matrix(rep(0,times=k*n),nrow=k)
    for(i in 1:n){
      for(c in 1:k){
        l[c,i]<-sum(dpois(y[,i],lambda=exp(coefs[,c]),log=TRUE))
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
    
    
    # print(pi) # print estimated cluster proportions
    
    
  }
  
  
  
  
  
  
  
  
  
  
  
  
  final_clusters<-rep(0,times=n)
  for(i in 1:n){
    final_clusters[i]<-which.max(finalwts[,i])
  }
  
  stats<-cluster.stats(d=d,true_clusters,final_clusters) #d is euclidean distance found earlier in cluster init.
  
  m<-rep(0,times=g)
  nondiscriminatory=rep(FALSE,times=g)
  true_nondiscriminatory=rep(FALSE,times=g)
  mean_across_clusters<-rowSums(coefs)/ncol(coefs)
  true_mean_across_clusters<-rowSums(b)/ncol(b)
  
  for(j in 1:g){
    if(all(abs(exp(b[j,])-exp(true_mean_across_clusters[j]))<7)){true_nondiscriminatory[j]=TRUE}     # threshold for nondiscriminatory gene: 1.5 diff from mean across clusters
    if(all(abs(exp(coefs[j,])-exp(mean_across_clusters[j]))<7)){nondiscriminatory[j]=TRUE}     # threshold for nondiscriminatory gene: 1.5 diff from mean across clusters
    for(c in 1:k){
      if(abs(exp(coefs[j,c])-exp(mean_across_clusters))>7){m[j]=m[j]+1} # nondiscriminatory threshold: away from mean by 7
    }
  }
  
  pred.nondiscriminatory<-mean(true_nondiscriminatory==nondiscriminatory)
  
  BIC=-2*Q[a]+log(n)*sum(m)
  
  
  
  result<-list(pi=pi,coefs=coefs,Q=Q[1:a],ARI=stats$corrected.rand,BIC=BIC,nondiscriminatory=pred.nondiscriminatory)
  #plot(Q[2:a])   # omit the first one due to instability
  return(result)
  
}
