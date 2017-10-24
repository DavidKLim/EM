library(caret)
library(mclust)
library(stats)
library(MASS)
library(fpc)
library(permute)
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
             size_factors=rep(1,times=ncol(y)) ,
             norm_y=y,
             true_clusters=NA){

n<-ncol(y)
g<-nrow(y)

# this makes it possible to have y=0 --> adds 0.1 to all y
y = y+0.1

# restructure data
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
#cls<-sample(1:k,n,replace=TRUE) #random initialization
  

########################## SIMULATION ONLY #############################
if(is.na(true_clusters)==FALSE){
  all_perms=allPerms(1:k)
  all_clusts=list()
  temp_clust<-rep(0,times=n)
  if(k==2){all_perms=matrix(all_perms,nrow=1)}
  
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
}
########################## SIMULATION ONLY #############################

  
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
  
coefs<-matrix(rep(0,times=g*k),nrow=g)
pi<-rep(0,times=k)
theta_list<-list()                # temporary to hold all K x K theta matrices
  
family=poisson(link="log")        # can specify family here
offset=log(size_factors)
#offset=rep(0,times=n)            # no offsets
  
IRLS_tol = 1E-7                   # Tolerance levels for embedded IRLS and Q fx in EM
maxit_IRLS=100

EM_tol = 1E-5
maxit_EM = 100
Q<-rep(0,times=maxit_EM)
  
lowerK<-0
  
  ########### M / E STEPS #########
  for(a in 1:maxit_EM){
    
  # M step
    
    dat[,"weights"]<-rep(as.vector(wts),times=g) # update weights column in dat
    # betaglm<-matrix(0,nrow=g,ncol=k)    # use to compare with glm procedure
    
    # IRWLS:
    for(j in 1:g){
      beta<-rep(0,times=k)

      if(a==1){
        for(c in 1:k){
          beta[c]<-log(mean(as.numeric(y[j,cls==c])))               # Initialize beta
        }
        
        theta<-matrix(rep(0,times=k^2),nrow=k)
        for(c in 1:k){
          for(cc in 1:k){
            if(theta[c,cc]>=tau){theta[c,cc]<-beta[c]-beta[cc]}             # Initialize theta
            else{theta[c,cc]<-soft_thresholding(beta[c]-beta[cc],lambda2)   # 
            }
          }
        }
      } else {
        beta<-coefs[j,]                                                   # Retrieve beta & theta from
        theta<-theta_list[[j]]                                            # previous iteration
      }
      
      temp<-matrix(rep(0,times=maxit_IRLS*k),nrow=maxit_IRLS)    # Temporarily store beta to test for convergence of IRLS
      dat_j<-dat[dat[,"g"]==j,]                                  # subset just the j'th gene
      
      for(i in 1:maxit_IRLS){
        
        temp[i,]<-beta
        eta <- 
          if(a==1 & i==1){
            matrix(rep(beta,times=n),nrow=n,byrow=TRUE)               # first initialization of eta
          }else if(a>1 & i==1){
            matrix(rep(beta,times=n),nrow=n,byrow=TRUE) + offset     # Retrieval of eta for IRLS (prev. beta + offset)
          }else{eta}
        
        for(c in 1:k){
          linkinv<-family$linkinv              # g^(-1) (eta) = mu
          mu.eta<-family$mu.eta         # g' = d(mu)/d(eta)
          variance<-family$variance
          
          mu = linkinv(eta)
          mu.eta.val = mu.eta(eta)
          
          dat_jc<-dat_j[dat_j[,"clusts"]==c,]    # subset j'th gene, c'th cluster
          
          good <- (dat_jc[,"weights"]>0) & (mu.eta.val[,c] != 0)
          
          trans_y <- (eta[,c] - offset)[good] + (dat_jc[,"count"][good] - mu[,c][good]) / mu.eta.val[,c][good]    # subtract size factor from transf. y
        
          w <- sqrt(dat_jc[,"weights"][good]*mu.eta.val[,c][good]^2/variance(mu[,c])[good])     # weights used in IRLS
          
          if(lambda1 != 0){            # Pan update
            beta[c]<-( (lambda1*((sum(beta)-beta[c]) + (sum(theta[c,])-theta[c,c])))  +  ((1/n)*sum(w*trans_y)) ) / ( (lambda1*(k-1)) + (1/n)*sum(w) )
          } else { beta[c]<-sum(w*trans_y)/sum(w) }
          
          #betaglm[j,c]<-log(glm(dat_jc[,"count"] ~ 1 + offset(offset), weights=dat_jc[,"weights"])$coef)   # glm update
          
          
          if(beta[c]<(-100)){
            warning(paste("Cluster",c,"Gene",j,"goes to -infinity"))
            beta[c] = -100
          }
          if(beta[c]>100){
            warning(paste("Cluster",c,"Gene",j,"goes to +infinity"))
            beta[c] = 100
          }
          
          eta[,c]<-beta[c] + offset      # add back size factors to eta
        }
        
        # update on theta
        for(c in 1:k){
          for(cc in 1:k){
            if(abs(theta[c,cc])>=tau){theta[c,cc]<-beta[c]-beta[cc]}                      # gTLP from Pan paper
            else{theta[c,cc]<-soft_thresholding(beta[c]-beta[cc],lambda2)}           # 
          }
        }
        
        # break conditions for IRLS
        if(i>1){
          if(sum((temp[i,]-temp[i-1,])^2)<IRLS_tol){
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
        l[c,i]<-sum(dpois(y[,i]-0.1,lambda=exp(coefs[,c]+offset[i]),log=TRUE))    # posterior log like, include size_factor of subj
      }  # correct for adding 0.1 earlier
    }
    
    # store and check Q function
    pt1<-(log(pi)%*%rowSums(wts))
    pt2<-sum(wts*l)
    Q[a]<-pt1+pt2
  
    # break condition for EM
    if(a>10){if(abs(Q[a]-Q[a-10])<EM_tol) {
      finalwts<-wts
      break
    }}
    
    
  # E step
    
    # update on weights
    logdenom = apply(log(pi) + l, 2,logsumexpc)
    for(c in 1:k){
      wts[c,]<-exp(log(pi[c])+l[c,]-logdenom)
    }
    
    if(a==maxit_EM){finalwts<-wts}
    # print(pi) # print estimated cluster proportions
    
    
    # if(any(rowSums(wts)==0)){
    #   print(paste("Empty cluster when K =",k,". Choose smaller K"))
    #   lowerK=1
    #   break
    # }
    for(i in 1:n){
      for(c in 1:k){
        if(wts[c,i]<1E-25){
          wts[c,i]=1E-25
        }
        if(wts[c,i]>(1-1E-25)){
          wts[c,i]=1E-25
        }
      }
    }
    
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
if(lowerK==1){BIC=.Machine$integer.max}      # set BIC as max (worst) if K too high
  
result<-list(pi=pi,
             coefs=coefs,
             Q=Q[1:a],
             BIC=BIC,
             nondiscriminatory=nondiscriminatory,
             final_clusters=final_clusters)
#plot(Q[2:a])   # omit the first one due to instability
return(result)
  
}
