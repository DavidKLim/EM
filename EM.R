library(mclust)
library(stats)
library(MASS)
library(fpc)

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

EM<-function(n,k,g,pi,b){

# simulate response data y matrix, g=20 by n=10000 matrix, with k=3 clusters #
y<-matrix(rep(0,times=g*n),nrow=g)      # initialize count matrix gxn #

z = rmultinom(n,1,pi)
for(j in 1:g){
  for(c in 1:k){
    y[j,z[c,]==1] = rpois(sum(z[c,]==1), lambda = exp(b[j,c]))
  }
}

# start EM

# initialization (HIERARCHICAL CLUSTERING)
#model<-Mclust(t(y),G=k)
#cls<-model$classification

d<-dist(t(y))
model<-hclust(d)
cls<-cutree(model,k=k)
#cls<-sample(1:k,n,replace=TRUE) #random initialization


wts<-matrix(rep(0,times=k*n),nrow=k)
for(c in 1:k){
  wts[c,]=(cls==c)^2
}
finalwts<-matrix(rep(0,times=k*n),nrow=k)


maxit = 100

coefs<-matrix(rep(0,times=g*k),nrow=g)
pi<-rep(0,times=k)
Q<-rep(0,times=maxit)

for(a in 1:maxit){
  
  # M step
  
  # estimate parameters by IRLS, gene by gene
  #for(j in 1:g){
  #  for(c in 1:k){
  #  fit = glm(y[j,] ~ 1, weights = wts[c,], family = poisson())   # automatically uses IRLS to find our maximization
  #  coefs[j,c]<-fit$coef
  #  }
  #}
  # IRWLS:
  lambda=.1
  for(j in 1:g){
    if(a>1) {eta<-coefs[j,]} else {eta<-rep(0.5,times=k)}
    temp<-matrix(rep(0,times=1000*k),nrow=1000)
    # w<-matrix(rep(0,times=n*k),nrow=k)    # not used b/c mu is directly coded in update
    for(i in 1:1000){
      mu=exp(eta)
      temp[i,]<-eta
      for (c in 1:k){
        nz<-eta[c]+wts[c,]*(y[j,]-mu[c])/mu[c]  # weighted observations (?) based on current weights w_jk
        # w[c,]<-rep(mu[c],times=length(nz)) # iterative weights for IRLS (not used)
        #m=0
        #for(d in 1:k){
        #  if(d != c){if(eta[c] != eta[d]) {m=m+abs(eta[c]-eta[d])}}
        #}
        eta[c]<-(mu[c]*mean(nz)+2*lambda*(sum(eta)-eta[c]))/(mu[c]+2*lambda*(k-1))       # gives negative value for last one (c=3)
        mu[c]<-exp(eta[c])
        temp[i,c]<-eta[c]
        if(i>1){
          if(sum((temp[i,]-temp[i-1,])^2)<1E-7){
            coefs[j,]<-eta
            break
          }
        }
      }
    }
  }
  
    
    
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

}

final_clusters<-rep(0,times=n)
true_clusters<-rep(0,times=n)
for(i in 1:n){
  final_clusters[i]<-which.max(finalwts[,i])
  true_clusters[i]<-which(z[,i]==1)
}

stats<-cluster.stats(d=d,true_clusters,final_clusters) #d is euclidean distance found earlier in cluster init.

result<-list(pi=pi,coefs=coefs,Q=Q[1:a],ARI=stats$corrected.rand)
plot(Q[2:a])   # omit the first one due to instability
return(result)
}



# parameters for simulation/algorithm
#n=20
#k=3
#g=200
#pi1=.3
#pi2=.3
#pi3=1-(pi1+pi2)
#pi<-c(pi1,pi2,pi3)

#b=matrix(rep(c(3,3.5,4),times=g),nrow=g,byrow=TRUE) # g x k matrix of betas
#b[1:50,]<-matrix(rep(c(2,3,2.5),times=50),nrow=50,byrow=TRUE)
#b[51:100,]<-matrix(rep(c(1,1.5,1.75),times=50),nrow=50,byrow=TRUE)
#b[101:150,]<-matrix(rep(c(1,1,3),times=50),nrow=50,byrow=TRUE)

# NEED TO ADD WAY TO INITIALIZE BETAS #

##### SIMULATIONS #####

#sim=100
#temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
#temp_coefs<-list()

#for(i in 1:sim){
#X<-EM(n=n,k=k,g=g,pi=pi,b=b)
#temp_pi[i,]<-X$pi
#temp_coefs[[i]]<-X$coefs
#}

#mean_pi<-colSums(t(apply(temp_pi,1,sort)))/sim

#for(i in 1:sim){
#  temp_coefs[[i]]<-t(apply(temp_coefs[[i]],1,sort))
#}
#mean_coefs<-Reduce('+',temp_coefs)/sim



#update on weights wrong
#Amax<-max(log(matrix(rep(pi,times=n),nrow=k))+l)
#logdenom=Amax+log(colSums(exp(log(matrix(rep(pi,times=n),nrow=k))+l-Amax)))
#for(c in 1:k){
#  wts[c,]<-exp(log(pi[c])+l[c,]-logdenom)
#}
