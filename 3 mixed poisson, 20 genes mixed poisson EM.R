library(mclust)

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
model<-Mclust(t(y),G=k)
cls<-model$classification

wts<-matrix(rep(0,times=k*n),nrow=k)
# note: unique(cls) should be 1:k, since k is specified
for(c in unique(cls)){
wts[c,]=(cls==c)^2
}

maxit = 200

coefs<-matrix(rep(0,times=g*k),nrow=g)
pi<-rep(0,times=k)
Q<-rep(0,times=maxit)
crit=1

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

# store and (ADD THIS: check stopping criterion)
pt1<-(log(pi)%*%rowSums(wts))
pt2<-sum(wts*l)
Q[a]<-pt1+pt2
if(a>10){if(abs(Q[a]-Q[a-10])<1E-5) break}

# E step

# update on weights [RECODE THIS]
Amax<-max(log(matrix(rep(pi,times=n),nrow=k))+l)
logdenom=Amax+log(colSums(exp(log(matrix(rep(pi,times=n),nrow=k))+l-Amax)))
wts[1,]<-exp(log(pi[1])+l[1,]-logdenom)
wts[2,]<-exp(log(pi[2])+l[2,]-logdenom)
wts[3,]<-exp(log(pi[3])+l[3,]-logdenom)

}

result<-list(pi=pi,coefs=coefs,Q=Q)
plot(Q[2:a])   # omit the first one due to instability
return(result)
}




# parameters for simulation/algorithm
n=1000
k=3
g=5
pi1=.3
pi2=.3
pi3=1-(pi1+pi2)
pi<-c(pi1,pi2,pi3)

b=matrix(rep(c(3,3.5,4),times=g),nrow=g,byrow=TRUE) # g x k matrix of betas
#b[1,]<-c(2,3,2.5)
#b[2,]<-c(1,1.5,1.75)
#b[3,]<-c(1,1,3)

# NEED TO ADD WAY TO INITIALIZE BETAS #

sim=100
temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
temp_coefs<-list()

for(i in 1:sim){
X<-EM(n=n,k=k,g=g,pi=pi,b=b)
temp_pi[i,]<-X$pi
temp_coefs[[i]]<-X$coefs
}

mean_pi<-colSums(t(apply(temp_pi,1,sort)))/sim

for(i in 1:sim){
  temp_coefs[[i]]<-t(apply(temp_coefs[[i]],1,sort))
}
mean_coefs<-Reduce('+',temp_coefs)/sim



# redo with pi=(.8,.1,.1)
pi=c(.8,.1,.1)
sim=100
temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
temp_coefs<-list()
for(i in 1:sim){
  X<-EM(n=n,k=k,g=g,pi=pi,b=b)
  temp_pi[i,]<-X$pi
  temp_coefs[[i]]<-X$coefs
}
mean_pi<-colSums(t(apply(temp_pi,1,sort)))/sim
for(i in 1:sim){
  temp_coefs[[i]]<-t(apply(temp_coefs[[i]],1,sort))
}
mean_coefs<-Reduce('+',temp_coefs)/sim



# redo with pi=(.3,.3,.4) and b1=(3,3.5,4) b2=(4,5,6)
pi=c(.3,.3,.4)
b[1,]=c(3,3.5,4)
b[2,]=c(4,5,6)
sim=100
temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
temp_coefs<-list()
for(i in 1:sim){
  X<-EM(n=n,k=k,g=g,pi=pi,b=b)
  temp_pi[i,]<-X$pi
  temp_coefs[[i]]<-X$coefs
}
mean_pi<-colSums(t(apply(temp_pi,1,sort)))/sim
for(i in 1:sim){
  temp_coefs[[i]]<-t(apply(temp_coefs[[i]],1,sort))
}
mean_coefs<-Reduce('+',temp_coefs)/sim





# redo with pi=(.8,.1,.1) and b1=(3,3.5,4) b2=(4,5,6)
pi=c(.8,.1,.1)
b[1,]=c(3,3.5,4)
b[2,]=c(4,5,6)
sim=100
temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
temp_coefs<-list()
for(i in 1:sim){
  X<-EM(n=n,k=k,g=g,pi=pi,b=b)
  temp_pi[i,]<-X$pi
  temp_coefs[[i]]<-X$coefs
}
mean_pi<-colSums(t(apply(temp_pi,1,sort)))/sim
for(i in 1:sim){
  temp_coefs[[i]]<-t(apply(temp_coefs[[i]],1,sort))
}
mean_coefs<-Reduce('+',temp_coefs)/sim





# redo with pi=(.8,.1,.1) and b1=(3,3.5,4) b2=(3.5,3.7,4)
pi=c(.8,.1,.1)
b[1,]=c(3,3.5,4)
b[2,]=c(3.5,3.7,4)
sim=100
temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
temp_coefs<-list()
for(i in 1:sim){
  X<-EM(n=n,k=k,g=g,pi=pi,b=b)
  temp_pi[i,]<-X$pi
  temp_coefs[[i]]<-X$coefs
}
mean_pi<-colSums(t(apply(temp_pi,1,sort)))/sim
for(i in 1:sim){
  temp_coefs[[i]]<-t(apply(temp_coefs[[i]],1,sort))
}
mean_coefs<-Reduce('+',temp_coefs)/sim




# redo with pi=(.8,.1,.1) and b1=(3,3.5,4) b2=(3.5,3.7,8)
pi=c(.8,.1,.1)
b[1,]=c(3,3.5,4)
b[2,]=c(3.5,3.7,8)
sim=100
temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
temp_coefs<-list()
for(i in 1:sim){
  X<-EM(n=n,k=k,g=g,pi=pi,b=b)
  temp_pi[i,]<-X$pi
  temp_coefs[[i]]<-X$coefs
}
mean_pi<-colSums(t(apply(temp_pi,1,sort)))/sim
for(i in 1:sim){
  temp_coefs[[i]]<-t(apply(temp_coefs[[i]],1,sort))
}
mean_coefs<-Reduce('+',temp_coefs)/sim


#################################################################
###################### failed attempts ##########################
#################################################################

# might cause underflow?

#for(c in 1:k){
#  denom=rep(0,n)
#  for(d in 1:k){
#    denom=denom+exp(l[d,])*pi[d]
#  }
#  wts[c,]=exp(l[c,])*pi[c]/denom
#}



# this is wrong (logsumexpc fx)
#logsumexpc=function(v){  
#  if(any(is.infinite(v))){ 
#    stop("infinite value in v\n") 
#  } 
#  if(length(v)==1){ return(v[1]) }  
#  sv = sort(v, decreasing=TRUE) 
#  res = sum(exp(sv[-1] - sv[1]))
#  lse = sv[1] + log(1+res)
#  lse 
#}
#pp1<-log(pi[1])+l[1,]
#pp2<-log(pi[2])+l[2,]
#pp3<-log(pi[3])+l[3,]
#probi1<-exp(pp1-logsumexpc(cbind(pp1,pp2,pp3)))
#probi2<-exp(pp2-logsumexpc(cbind(pp1,pp2,pp3)))
#probi3<-1-(probi1+probi2)


# logsumexp trick online

#Amax<-rep(max(log(pi*exp(l))),n)
#sumexp=rep(0,n)
#for(d in 1:k){
#  sumexp=sumexp+exp(log(pi[d]*exp(l[d,]))-Amax)
#}
#logdenom=Amax+log(sumexp)
#for(c in 1:k){
#  wts[c,] = exp(log(pi[c]*exp(l[c,])) - logdenom)
#}
