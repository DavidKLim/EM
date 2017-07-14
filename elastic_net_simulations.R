setwd("/netscr/deelim")
#setwd("C:/Users/David/Desktop/Research/GitHub/EM")
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
# alpha_search=seq(0,1,by=0.25) # 5 elts
alpha_search=1    # rerun with alpha=1
lambda1_search=c(0.01,0.05, 0.1, 0.2, 0.5, 0.8, 1) #
lambda2_search=c(0.01,0.05, 0.1, 0.2, 0.5, 0.8, 1) #

list_BIC=matrix(0,nrow=length(alpha_search)*length(lambda1_search)*length(lambda2_search),ncol=4) #matrix of BIC's: alpha x lambda

list_BIC[,1]=rep(alpha_search,each=length(lambda1_search)*length(lambda2_search))
list_BIC[,2]=rep(rep(lambda1_search,each=length(lambda2_search)),times=length(alpha_search))
list_BIC[,3]=rep(lambda2_search,times=length(alpha_search)*length(lambda1_search))

for(aa in 1:nrow(list_BIC)){
  set.seed(aa)
  list_BIC[aa,4]<-EM(n=n,k=k,g=g,init_pi=pi,b=b,alpha=list_BIC[aa,1],lambda1=list_BIC[aa,2],lambda2=list_BIC[aa,3])$BIC
  print(list_BIC[aa,])
}

max_index<-which(list_BIC[,4]==min(list_BIC[,4]))

max_alpha<-list_BIC[max_index,1]
max_lambda1<-list_BIC[max_index,2]
max_lambda2<-list_BIC[max_index,3]

num_nondiscr<-rep(0,times=sim)

cl_accuracy<-rep(0,times=sim)
sens<-rep(0,times=sim)
false_pos<-rep(0,times=sim)


#simulations
for(i in 1:sim){
  set.seed(i)
  X<-EM(n=n,k=k,g=g,init_pi=pi,b=b,alpha=max_alpha,lambda1=max_lambda1,lambda2=max_lambda2)
  temp_pi[i,]<-X$pi
  temp_coefs[[i]]<-X$coefs
  ARI[i]<-X$ARI
  num_nondiscr[i]<-mean(X$nondiscriminatory)
  cl_accuracy[i]<-X$cl_accuracy
  sens[i]<-X$sens
  false_pos[i]<-X$false_pos
}

mean_pi<-colSums(temp_pi)/sim
mean_coefs<-Reduce('+',temp_coefs)/sim
coef1_means<-colSums(mean_coefs)/nrow(mean_coefs)

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
print(paste("Mean Cluster Accuracy =",mean(cl_accuracy)))
print(paste("Mean sensitivity =", mean(sens)))
print(paste("Mean false positive rate =", mean(false_pos)))
sink()






















# CASE 2: 50% nondiscriminatory
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
# alpha_search=seq(0,1,by=0.25) # 5 elts
alpha_search=1    # rerun with alpha=1
lambda1_search=c(0.01,0.05, 0.1, 0.2, 0.5, 0.8, 1) #
lambda2_search=c(0.01,0.05, 0.1, 0.2, 0.5, 0.8, 1) #

list_BIC=matrix(0,nrow=length(alpha_search)*length(lambda1_search)*length(lambda2_search),ncol=4) #matrix of BIC's: alpha x lambda

list_BIC[,1]=rep(alpha_search,each=length(lambda1_search)*length(lambda2_search))
list_BIC[,2]=rep(rep(lambda1_search,each=length(lambda2_search)),times=length(alpha_search))
list_BIC[,3]=rep(lambda2_search,times=length(alpha_search)*length(lambda1_search))

for(aa in 1:nrow(list_BIC)){
  set.seed(aa)
  list_BIC[aa,4]<-EM(n=n,k=k,g=g,init_pi=pi,b=b,alpha=list_BIC[aa,1],lambda1=list_BIC[aa,2],lambda2=list_BIC[aa,3])$BIC
  print(list_BIC[aa,])
}

max_index<-which(list_BIC[,4]==min(list_BIC[,4]))

max_alpha<-list_BIC[max_index,1]
max_lambda1<-list_BIC[max_index,2]
max_lambda2<-list_BIC[max_index,3]

num_nondiscr<-rep(0,times=sim)

cl_accuracy<-rep(0,times=sim)
sens<-rep(0,times=sim)
false_pos<-rep(0,times=sim)

#simulations
for(i in 1:sim){
  set.seed(i)
  X<-EM(n=n,k=k,g=g,init_pi=pi,b=b,alpha=max_alpha,lambda1=max_lambda1,lambda2=max_lambda2)
  temp_pi[i,]<-X$pi
  temp_coefs[[i]]<-X$coefs
  ARI[i]<-X$ARI
  num_nondiscr[i]<-mean(X$nondiscriminatory)
  cl_accuracy[i]<-X$cl_accuracy
  sens[i]<-X$sens
  false_pos[i]<-X$false_pos
}

mean_pi<-colSums(temp_pi)/sim
mean_coefs<-Reduce('+',temp_coefs)/sim
coef1_means<-colSums(mean_coefs[1:100,])/100
coef2_means<-colSums(mean_coefs[101:200,])/100

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
print(paste("Mean Cluster Accuracy =",mean(cl_accuracy)))
print(paste("Mean sensitivity =", mean(sens)))
print(paste("Mean false positive rate =", mean(false_pos)))
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
# alpha_search=seq(0,1,by=0.25) # 5 elts
alpha_search=1    # rerun with alpha=1
lambda1_search=c(0.01,0.05, 0.1, 0.2, 0.5, 0.8, 1) #
lambda2_search=c(0.01,0.05, 0.1, 0.2, 0.5, 0.8, 1) #

list_BIC=matrix(0,nrow=length(alpha_search)*length(lambda1_search)*length(lambda2_search),ncol=4) #matrix of BIC's: alpha x lambda

list_BIC[,1]=rep(alpha_search,each=length(lambda1_search)*length(lambda2_search))
list_BIC[,2]=rep(rep(lambda1_search,each=length(lambda2_search)),times=length(alpha_search))
list_BIC[,3]=rep(lambda2_search,times=length(alpha_search)*length(lambda1_search))

for(aa in 1:nrow(list_BIC)){
  set.seed(aa)
  list_BIC[aa,4]<-EM(n=n,k=k,g=g,init_pi=pi,b=b,alpha=list_BIC[aa,1],lambda1=list_BIC[aa,2],lambda2=list_BIC[aa,3])$BIC
  print(list_BIC[aa,])
}

max_index<-which(list_BIC[,4]==min(list_BIC[,4]))

max_alpha<-list_BIC[max_index,1]
max_lambda1<-list_BIC[max_index,2]
max_lambda2<-list_BIC[max_index,3]

num_nondiscr<-rep(0,times=sim)

cl_accuracy<-rep(0,times=sim)
sens<-rep(0,times=sim)
false_pos<-rep(0,times=sim)

#simulations

for(i in 1:sim){
  set.seed(i)
  X<-EM(n=n,k=k,g=g,init_pi=pi,b=b,alpha=max_alpha,lambda1=max_lambda1,lambda2=max_lambda2)
  temp_pi[i,]<-X$pi
  temp_coefs[[i]]<-X$coefs
  ARI[i]<-X$ARI
  num_nondiscr[i]<-mean(X$nondiscriminatory)
  cl_accuracy[i]<-X$cl_accuracy
  sens[i]<-X$sens
  false_pos[i]<-X$false_pos
}

mean_pi<-colSums(temp_pi)/sim
mean_coefs<-Reduce('+',temp_coefs)/sim
coef1_means<-colSums(mean_coefs[1:20,])/20
coef2_means<-colSums(mean_coefs[21:200,])/180

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
print(paste("Mean Cluster Accuracy =",mean(cl_accuracy)))
print(paste("Mean sensitivity =", mean(sens)))
print(paste("Mean false positive rate =", mean(false_pos)))
sink()