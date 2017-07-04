setwd("/netscr/deelim")
#setwd("C:/Users/David/Desktop/Research/Coding EM")
source("EM.R")
library(MASS)

n=20
k=3
g=200
pi=c(.5,.3,.2)
sigma=diag(k)
b=matrix(rep(0,times=k*g),nrow=g,byrow=TRUE) # initialize betas

#### SIMULATIONS ####

# CASE 1: distinct means, uniform across genes
b[1:200,]<-mvrnorm(200,mu=c(3,3.5,4),sigma) # betas with very distinct means

sim=100
temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
temp_coefs<-list()
ARI<-rep(0,times=sim)

#grid search
alpha_search=seq(0,1,by=0.25) # 5 elts
lambda1_search=seq(0.1,1,by=0.1) #lambda1 can't = 0. 10 elts
lambda2_search=c(0.01,0.05,0.1,0.2,1) # 5 elts

list_BIC=matrix(0,nrow=length(alpha_search)*length(lambda1_search)*length(lambda2_search),ncol=4) #matrix of BIC's: alpha x lambda

list_BIC[,1]=rep(alpha_search,each=length(lambda1_search)*length(lambda2_search))
list_BIC[,2]=rep(rep(lambda1_search,each=length(lambda2_search)),times=length(alpha_search))
list_BIC[,3]=rep(lambda2_search,times=length(alpha_search)*length(lambda1_search))

for(aa in 1:nrow(list_BIC)){
  list_BIC[aa,4]<-EM(n=n,k=k,g=g,init_pi=pi,b=b,alpha=list_BIC[aa,1],lambda1=list_BIC[aa,2],lambda2=list_BIC[aa,3])$BIC
  print(list_BIC[aa,])
}

max_index<-which(list_BIC[,4]==min(list_BIC[,4]))

max_alpha<-list_BIC[max_index,1]
max_lambda1<-list_BIC[max_index,2]
max_lambda2<-list_BIC[max_index,3]

num_nondiscr<-rep(0,times=sim)


#simulations
for(i in 1:sim){
  X<-EM(n=n,k=k,g=g,init_pi=pi,b=b,alpha=max_alpha,lambda1=max_lambda1,lambda2=max_lambda2)
  temp_pi[i,]<-X$pi
  temp_coefs[[i]]<-X$coefs
  ARI[i]<-X$ARI
  num_nondiscr[i]<-X$nondiscriminatory
}

mean_pi<-colSums(t(apply(temp_pi,1,sort)))/sim # sort in increasing order, then mean #
mean_pi<-mean_pi[order(mean_pi)][rank(pi)] # order of true pi

# array of gxk matrices (100 (numsims) of them) of coefficients
sorted_coefs<-list()
for(i in 1:sim){
  sorted_coefs[[i]]<-temp_coefs[[i]][,order(temp_pi[i,])] # sort by increasing order of pi (to align --> take mean) #
}   # trickier when varying each gene mean

mean_coefs<-Reduce('+',sorted_coefs)/sim # gxk matrix of mean coef's across simulations #
sorted_means<-matrix(rep(0,times=g*k),nrow=g)
for(j in 1:g){
  sorted_means[j,]<-mean_coefs[j,rank(pi)] # order each row according to order of true pi #
}
coef1_means<-colSums(sorted_means)/nrow(sorted_means)

mean_ARI<-mean(ARI)

SSE_estims<-sum((pi-mean_pi)^2)+sum((b-sorted_means)^2)




# output #
sink("elastic net.txt")
print("CASE 1: b ~ MVN(mu=(3, 3.5, 4)) all genes: 0% nondiscriminatory")
print(paste("n =",n,", k =",k,", g =",g,", pi =",pi[1],pi[2],pi[3]))
print(paste("simulated cluster means: MVN(mu=(3, 3.5, 4)"))
print(paste("mean_pi =",mean_pi[1],mean_pi[2],mean_pi[3])) # mean of estimated pi
print(paste("mean_coefs =",coef1_means[1],coef1_means[2],coef1_means[3])) # in varying means, needs work
print(paste("mean_ARI =",mean_ARI)) # mean of corrected rand index
print(paste("SSE_estimates =",SSE_estims)) # sum of square errors of all estimated parameters (pi and coefs)
print(paste("mean % of nondiscriminatory genes =",mean(num_nondiscr)))
print(paste("final (alpha,lambda1,lambda2) =",max_alpha,max_lambda1,max_lambda2))
sink()






















# CASE 3: 50% nondiscriminatory
n=20
k=3
g=200
pi=c(.5,.3,.2)
sigma=diag(k)

b=matrix(rep(c(1,1,1),times=g),nrow=g,byrow=TRUE) # initialize
b[1:100,]<-mvrnorm(100,mu=c(1,2,3),sigma)
b[101:200,]<-mvrnorm(100,mu=c(2,2,2),sigma) # 50% nondiscriminatory

##### SIMULATIONS #####
# 100% nondiscriminatory #
sim=100
temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
temp_coefs<-list()
ARI<-rep(0,times=sim)


#grid search
alpha_search=seq(0,1,by=0.25) # 5 elts
lambda1_search=seq(0.1,1,by=0.1) #lambda1 can't = 0. 10 elts
lambda2_search=c(0.01,0.05,0.1,0.2,1) # 5 elts

list_BIC=matrix(0,nrow=length(alpha_search)*length(lambda1_search)*length(lambda2_search),ncol=4) #matrix of BIC's: alpha x lambda

list_BIC[,1]=rep(alpha_search,each=length(lambda1_search)*length(lambda2_search))
list_BIC[,2]=rep(rep(lambda1_search,each=length(lambda2_search)),times=length(alpha_search))
list_BIC[,3]=rep(lambda2_search,times=length(alpha_search)*length(lambda1_search))

for(aa in 1:nrow(list_BIC)){
  list_BIC[aa,4]<-EM(n=n,k=k,g=g,init_pi=pi,b=b,alpha=list_BIC[aa,1],lambda1=list_BIC[aa,2],lambda2=list_BIC[aa,3])$BIC
  print(list_BIC[aa,])
}

max_index<-which(list_BIC[,4]==min(list_BIC[,4]))

max_alpha<-list_BIC[max_index,1]
max_lambda1<-list_BIC[max_index,2]
max_lambda2<-list_BIC[max_index,3]

num_nondiscr<-rep(0,times=sim)

#simulations
for(i in 1:sim){
  X<-EM(n=n,k=k,g=g,init_pi=pi,b=b,alpha=max_alpha,lambda1=max_lambda1,lambda2=max_lambda2)
  temp_pi[i,]<-X$pi
  temp_coefs[[i]]<-X$coefs
  ARI[i]<-X$ARI
  num_nondiscr[i]<-X$nondiscriminatory
}

mean_pi<-colSums(t(apply(temp_pi,1,sort)))/sim # sort in increasing order, then mean #
mean_pi<-mean_pi[order(mean_pi)][rank(pi)] # order of true pi

# array of gxk matrices (100 (numsims) of them) of coefficients
sorted_coefs<-list()
for(i in 1:sim){
  sorted_coefs[[i]]<-temp_coefs[[i]][,order(temp_pi[i,])] # sort by increasing order of pi (to align --> take mean) #
}   # trickier when varying each gene mean

mean_coefs<-Reduce('+',sorted_coefs)/sim # gxk matrix of mean coef's across simulations #
sorted_means<-matrix(rep(0,times=g*k),nrow=g)
for(j in 1:g){
  sorted_means[j,]<-mean_coefs[j,rank(pi)] # order each row according to order of true pi #
}
coef1_means<-colSums(sorted_means[1:100,])/nrow(sorted_means[1:100,])
coef2_means<-colSums(sorted_means[101:200,])/nrow(sorted_means[101:200,])

mean_ARI<-mean(ARI)

SSE_estims<-sum((pi-mean_pi)^2)+sum((b-sorted_means)^2)

# output #
sink("elastic net.txt",append=TRUE)
print("CASE 2: 50% nondiscriminatory")
print(paste("n =",n,", k =",k,", g =",g,", pi =",pi[1],pi[2],pi[3]))
print(paste("simulated cluster means: first 100: MVN(mu=(1,2,3),sigma=I), next 100: MVN(mu=(2,2,2),sigma=I)"))
print(paste("mean_pi =",mean_pi[1],mean_pi[2],mean_pi[3])) # mean of estimated pi
print(paste("mean_coefs1 =",coef1_means[1],coef1_means[2],coef1_means[3])) # in varying means, needs work
print(paste("mean_coefs2 =",coef2_means[1],coef2_means[2],coef2_means[3]))
print(paste("mean_ARI =",mean_ARI)) # mean of corrected rand index
print(paste("SSE_estimates =",SSE_estims)) # sum of square errors of all estimated parameters (pi and coefs)
print(paste("mean % of nondiscriminatory genes =",mean(num_nondiscr)))
print(paste("final (alpha,lambda1,lambda2) =",max_alpha,max_lambda1,max_lambda2))
sink()
















# CASE 3: 90% nondiscriminatory
n=20
k=3
g=200
pi=c(.5,.3,.2)
sigma=diag(k)

b=matrix(rep(c(1,1,1),times=g),nrow=g,byrow=TRUE) # initialize
b[1:20,]<-mvrnorm(20,mu=c(1,2,3),sigma)
b[21:200,]<-mvrnorm(180,mu=c(2,2,2),sigma) # 90% nondiscriminatory

##### SIMULATIONS #####
sim=100
temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
temp_coefs<-list()
ARI<-rep(0,times=sim)


#grid search
alpha_search=seq(0,1,by=0.25) # 5 elts
lambda1_search=seq(0.1,1,by=0.1) #lambda1 can't = 0. 10 elts
lambda2_search=c(0.01,0.05,0.1,0.2,1) # 5 elts

list_BIC=matrix(0,nrow=length(alpha_search)*length(lambda1_search)*length(lambda2_search),ncol=4) #matrix of BIC's: alpha x lambda

list_BIC[,1]=rep(alpha_search,each=length(lambda1_search)*length(lambda2_search))
list_BIC[,2]=rep(rep(lambda1_search,each=length(lambda2_search)),times=length(alpha_search))
list_BIC[,3]=rep(lambda2_search,times=length(alpha_search)*length(lambda1_search))

for(aa in 1:nrow(list_BIC)){
  list_BIC[aa,4]<-EM(n=n,k=k,g=g,init_pi=pi,b=b,alpha=list_BIC[aa,1],lambda1=list_BIC[aa,2],lambda2=list_BIC[aa,3])$BIC
  print(list_BIC[aa,])
}

max_index<-which(list_BIC[,4]==min(list_BIC[,4]))

max_alpha<-list_BIC[max_index,1]
max_lambda1<-list_BIC[max_index,2]
max_lambda2<-list_BIC[max_index,3]

num_nondiscr<-rep(0,times=sim)

#simulations

for(i in 1:sim){
  X<-EM(n=n,k=k,g=g,init_pi=pi,b=b,alpha=max_alpha,lambda1=max_lambda1,lambda2=max_lambda2)
  temp_pi[i,]<-X$pi
  temp_coefs[[i]]<-X$coefs
  ARI[i]<-X$ARI
  num_nondiscr[i]<-X$nondiscriminatory
}

mean_pi<-colSums(t(apply(temp_pi,1,sort)))/sim # sort in increasing order, then mean #
mean_pi<-mean_pi[order(mean_pi)][rank(pi)] # order of true pi

# array of gxk matrices (100 (numsims) of them) of coefficients
sorted_coefs<-list()
for(i in 1:sim){
  sorted_coefs[[i]]<-temp_coefs[[i]][,order(temp_pi[i,])] # sort by increasing order of pi (to align --> take mean) #
}   # trickier when varying each gene mean

mean_coefs<-Reduce('+',sorted_coefs)/sim # gxk matrix of mean coef's across simulations #
sorted_means<-matrix(rep(0,times=g*k),nrow=g)
for(j in 1:g){
  sorted_means[j,]<-mean_coefs[j,rank(pi)] # order each row according to order of true pi #
}
coef1_means<-colSums(sorted_means[1:20,])/nrow(sorted_means[1:20,])
coef2_means<-colSums(sorted_means[21:200,])/nrow(sorted_means[21:200,])

mean_ARI<-mean(ARI)

SSE_estims<-sum((pi-mean_pi)^2)+sum((b-sorted_means)^2)

# output #
sink("elastic net.txt",append=TRUE)
print("CASE 3: 90% nondiscriminatory")
print(paste("n =",n,", k =",k,", g =",g,", pi =",pi[1],pi[2],pi[3]))
print(paste("simulated cluster means: first 20: MVN(mu=(1,2,3),sigma=I), next 180: MVN(mu=(2,2,2),sigma=I)"))
print(paste("mean_pi =",mean_pi[1],mean_pi[2],mean_pi[3])) # mean of estimated pi
print(paste("mean_coefs1 =",coef1_means[1],coef1_means[2],coef1_means[3])) # in varying means, needs work
print(paste("mean_coefs2 =",coef2_means[1],coef2_means[2],coef2_means[3]))
print(paste("mean_ARI =",mean_ARI)) # mean of corrected rand index
print(paste("SSE_estimates =",SSE_estims)) # sum of square errors of all estimated parameters (pi and coefs)
print(paste("mean % of nondiscriminatory genes =",mean(num_nondiscr)))
print(paste("final (alpha,lambda1,lambda2) =",max_alpha,max_lambda1,max_lambda2))
sink()