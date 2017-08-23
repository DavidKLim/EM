setwd("/netscr/deelim")
#setwd("C:/Users/David/Desktop/Research/GitHub/EM")
source("non-simulation EM.R")   # Pan method
library(MASS)

n=20
k=4
g=500
pi=c(0.3,0.3,0.2,0.2)
sigma=diag(k)
b=matrix(rep(0,times=k*g),nrow=g,byrow=TRUE) # initialize betas


#### SIMULATIONS ####

# CASE 1: distinct means, uniform across genes
b[1:500,]<-mvrnorm(500,mu=c(10,10.5,11,11.5),sigma) # betas with very distinct means
b[1:250,]<-mvrnorm(250,mu=c(9.5,9.5,9.5,9.5),sigma)

b[1:500,]<-matrix(rep(c(10,10.5,11,11.5),times=500),nrow=500,byrow=TRUE) # Fixing the means to ensure no nondiscriminatory cases
b[1:250,]<-matrix(rep(c(9,9,9,9),times=250),nrow=250)

sim_size_factors<-rep(1,times=n) ### Size Factors
sim_size_factors<-seq(from=0.75, to=2, length.out=n)

simulate_data=function(n,k,g,init_pi,b){
  y<-matrix(rep(0,times=g*n),nrow=g)      # initialize count matrix gxn #
  # Prepare new flattened data
  z = rmultinom(n,1,init_pi)
  # while(any(rowSums(z)==0)){z=rmultinom(n,1,init_pi)}   # makes sure that no one cluster simulated @ 0 membership (only good for simulations)
  for(j in 1:g){
    for(c in 1:k){
      y[j,z[c,]==1] = rpois(sum(z[c,]==1), lambda = sim_size_factors[z[c,]==1]*exp(b[j,c]))
    }
  }
  result<-list(y=y,z=z)
  return(result)
}
sim.dat<-simulate_data(n=n,k=k,g=g,init_pi=pi,b=b)
y<-sim.dat$y
z<-sim.dat$z

true_clusters<-rep(0,times=n)
for(i in 1:n){
  true_clusters[i]<-which(z[,i]==1)
}



library("DESeq2")
library("genefilter")
row_names<-paste("gene",seq(g))
col_names<-paste("subj",seq(n))

cts<-as.matrix(y)
rownames(cts)<-row_names
colnames(cts)<-col_names
coldata<-matrix(paste("cl",true_clusters,sep=""),nrow=n)

rownames(coldata)<-colnames(cts)
colnames(coldata)<-"cluster"

dds<-DESeqDataSetFromMatrix(countData = cts,
                            colData = coldata,
                            design = ~ 1)
DESeq_dds<-DESeq(dds)
size_factors<-estimateSizeFactors(dds)$sizeFactor


norm_y<-counts(DESeq_dds,normalized=TRUE)



######### Order Selection (using unpenalized model) ##########
source("C:/Users/David/Desktop/Research/GitHub/EM/unpenalized Pan EM.R")
K_search=c(2:15)
list_BIC=matrix(0,nrow=length(K_search),ncol=2)
list_BIC[,1]=K_search
for(aa in 1:nrow(list_BIC)){
  list_BIC[aa,2]<-EM(y=y,k=list_BIC[aa,1],size_factors=size_factors)$BIC
  print(list_BIC[aa,])
}

max_k=list_BIC[which(list_BIC[,2]==min(list_BIC[,2])),1]





########## PAN ##########
source("C:/Users/David/Desktop/Research/GitHub/EM/Pan EM.R")
lambda1_search=seq(from=0.1,to=20,length.out=20)
lambda2_search=seq(from=0.1,to=20,length.out=20)
tau_search=seq(from=1,to=2,length.out=10) # nullifies tau param

list_BIC=matrix(0,nrow=length(lambda1_search)*length(lambda2_search)*length(tau_search),ncol=4) #matrix of BIC's: lambda1 and lambda2 and K, 49*5 combinations

list_BIC[,1]=rep(lambda1_search,each=length(lambda2_search)*length(tau_search))
list_BIC[,2]=rep(rep(lambda2_search,each=length(tau_search)),times=length(lambda1_search))
list_BIC[,3]=rep(tau_search,times=length(lambda1_search)*length(lambda2_search))

for(aa in 1:nrow(list_BIC)){
  list_BIC[aa,4]<-EM(y=y,k=max_k,tau=list_BIC[aa,3],lambda1=list_BIC[aa,1],lambda2=list_BIC[aa,2],size_factors=size_factors)$BIC
  print(list_BIC[aa,])
}

max_index<-which(list_BIC[,4]==min(list_BIC[,4]))

max_tau<-list_BIC[max_index,3]
max_lambda1<-list_BIC[max_index,1]
max_lambda2<-list_BIC[max_index,2]







######### GLMNET ########
source("C:/Users/David/Desktop/Research/GitHub/EM/group lasso EM.R")
alpha_search=seq(from=0,to=1,by=0.2)
lambda_search=seq(from=0,to=5,by=0.5)
list_BIC=matrix(0,nrow=length(alpha_search)*length(lambda_search),ncol=3)
list_BIC[,1]=rep(alpha_search,each=length(lambda_search))
list_BIC[,2]=rep(lambda_search,times=length(alpha_search))

for(aa in 1:nrow(list_BIC)){
  list_BIC[aa,3]<-EM(y=y,k=max_k,alpha=list_BIC[aa,1],lambda=list_BIC[aa,2],size_factors<-size_factors)$BIC
  print(list_BIC[aa,])
}

max_index<-which(list_BIC[,3]==min(list_BIC[,3]))

max_alpha<-list_BIC[max_index,1]
max_lambda<-list_BIC[max_index,2]













sim=100
temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
temp_coefs<-list()
ARI<-rep(0,times=sim)
num_nondiscr<-rep(0,times=sim)
cl_accuracy<-rep(0,times=sim)
sens<-rep(0,times=sim)
false_pos<-rep(0,times=sim)


#simulations
for(i in 1:sim){
  set.seed(i)
  X<-EM(n=n,k=k,g=g,init_pi=pi,b=b,lambda1=max_lambda1,lambda2=max_lambda2)
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
print(paste("final (lambda1,lambda2) =",max_lambda1,max_lambda2))
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
  list_BIC[aa,4]<-EM(n=n,k=k,g=g,init_pi=pi,b=b,lambda1=list_BIC[aa,2],lambda2=list_BIC[aa,3])$BIC
  print(list_BIC[aa,])
}

max_index<-which(list_BIC[,4]==min(list_BIC[,4]))

max_lambda1<-list_BIC[max_index,2]
max_lambda2<-list_BIC[max_index,3]

num_nondiscr<-rep(0,times=sim)

cl_accuracy<-rep(0,times=sim)
sens<-rep(0,times=sim)
false_pos<-rep(0,times=sim)

#simulations
for(i in 1:sim){
  set.seed(i)
  X<-EM(n=n,k=k,g=g,init_pi=pi,b=b,lambda1=max_lambda1,lambda2=max_lambda2)
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
print(paste("final (lambda1,lambda2) =",max_lambda1,max_lambda2))
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
  list_BIC[aa,4]<-EM(n=n,k=k,g=g,init_pi=pi,b=b,lambda1=list_BIC[aa,2],lambda2=list_BIC[aa,3])$BIC
  print(list_BIC[aa,])
}

max_index<-which(list_BIC[,4]==min(list_BIC[,4]))

max_lambda1<-list_BIC[max_index,2]
max_lambda2<-list_BIC[max_index,3]

num_nondiscr<-rep(0,times=sim)

cl_accuracy<-rep(0,times=sim)
sens<-rep(0,times=sim)
false_pos<-rep(0,times=sim)

#simulations

for(i in 1:sim){
  set.seed(i)
  X<-EM(n=n,k=k,g=g,init_pi=pi,b=b,lambda1=max_lambda1,lambda2=max_lambda2)
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
print(paste("final (lambda1,lambda2) =",max_lambda1,max_lambda2))
print(paste("Mean Cluster Accuracy =",mean(cl_accuracy)))
print(paste("Mean sensitivity =", mean(sens)))
print(paste("Mean false positive rate =", mean(false_pos)))
sink()