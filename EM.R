library(mclust)
library(stats)


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


wts<-matrix(rep(0,times=k*n),nrow=k)
for(c in 1:k){
  wts[c,]=(cls==c)^2
}


maxit = 200

coefs<-matrix(rep(0,times=g*k),nrow=g)
pi<-rep(0,times=k)
Q<-rep(0,times=maxit)

for(a in 1:maxit){

# M step

# estimate parameters by IRLS, gene by gene
for(j in 1:g){
  for(c in 1:k){
  fit = glm(y[j,] ~ 1, weights = wts[c,], family = poisson())   # automatically uses IRLS to find our maximization
  coefs[j,c]<-fit$coef
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
if(a>10){if(abs(Q[a]-Q[a-10])<1E-5) break}

# E step

# update on weights
Amax<-max(log(matrix(rep(pi,times=n),nrow=k))+l)
logdenom=Amax+log(colSums(exp(log(matrix(rep(pi,times=n),nrow=k))+l-Amax)))
for(c in 1:k){
  wts[c,]<-exp(log(pi[c])+l[c,]-logdenom)
}

}

result<-list(pi=pi,coefs=coefs,Q=Q)
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

