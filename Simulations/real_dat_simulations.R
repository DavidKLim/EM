#setwd("/netscr/deelim")
setwd("C:/Users/David/Desktop/Research/EM")
source("non-simulation EM.R")   # Pan method
library(MASS)



# Simulations to choose K
sim=100
choose_k<-rep(0,times=sim)
n=20
g=100
k=4
pi=c(0.2,0.4,0.3,0.1)
sigma=diag(k)
b=matrix(rep(0,times=k*g),nrow=g,byrow=TRUE) # initialize betas
b[1:100,]<-matrix(rep(c(10,10.5,11,9.5),times=100),nrow=100,byrow=TRUE) # Fixing the means to ensure no nondiscriminatory cases
b[1:50,]<-matrix(rep(c(9.5,9.5,9.5,9.5),times=50),nrow=50)

for(ii in 1:sim){
  simulate_data=function(n,k,g,init_pi,b){
    y<-matrix(rep(0,times=g*n),nrow=g)      # initialize count matrix gxn #
    # Prepare new flattened data
    z = rmultinom(n,1,init_pi)
    # while(any(rowSums(z)==0)){z=rmultinom(n,1,init_pi)}   # makes sure that no one cluster simulated @ 0 membership (only good for simulations)
    for(j in 1:g){
      for(c in 1:k){
        y[j,z[c,]==1] = rpois(sum(z[c,]==1), lambda = exp(b[j,c])*true_size_factors)
      }
    }
    result<-list(y=y,z=z)
    return(result)
  }
  sim.dat<-simulate_data(n=n,k=k,g=g,init_pi=pi,b=b)
  y<-sim.dat$y+1
  z<-sim.dat$z
  
  true_clusters<-rep(0,times=n)
  for(i in 1:n){
    true_clusters[i]<-which(z[,i]==1)
  }
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
   # scaled_y<-y
   # for(i in 1:n){
   #   scaled_y[,i]<-y[,i]/size_factors[i]
   # }
  
  ######### Order Selection (using unpenalized model) ##########
  source("C:/Users/David/Desktop/Research/EM/Pan EM.R")
  #source("C:/Users/David/Desktop/Research/EM/unpenalized EM.R")
  K_search=c(2:8)
  list_BIC=matrix(0,nrow=length(K_search),ncol=2)
  list_BIC[,1]=K_search
  
  print(paste("Iteration",ii,":"))
  for(aa in 1:nrow(list_BIC)){
     #nam <- paste("Xpen",list_BIC[aa,1],sep="")
     #assign(nam,EM(y=y,k=list_BIC[aa,1],lambda1=0,lambda2=0,tau=0,size_factors=size_factors))
    list_BIC[aa,2]<-EM(y=y,k=list_BIC[aa,1],lambda1=0,lambda2=0,tau=0,size_factors=size_factors)$BIC   # no penalty Pan
    #list_BIC[aa,2]<-EM(y=y,k=list_BIC[aa,1],size_factors=size_factors)$BIC       # unpenalized (not Pan)
    print(list_BIC[aa,])
  }
  
  max_k=list_BIC[which(list_BIC[,2]==min(list_BIC[,2])),1]
  choose_k[ii]<-max_k
}

table(choose_k)





# library("optCluster")
# opt.cl<-optCluster(round(scaled_y,0),2:8,clMethods="em.poisson",countData=TRUE)









########## PAN ##########
source("C:/Users/David/Desktop/Research/EM/Pan EM.R")
lambda1_search=seq(from=0.1,to=2,length.out=10)
lambda2_search=seq(from=0.1,to=2,length.out=10)
tau_search=seq(from=1,to=2,length.out=5) # nullifies tau param

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
# 
# 
# 
# 
# 
# 
# 
# ######### GLMNET ########
# source("C:/Users/David/Desktop/Research/EM/group lasso EM.R")
# alpha_search=seq(from=0,to=1,by=0.2)
# lambda_search=seq(from=0,to=5,by=0.5)
# list_BIC=matrix(0,nrow=length(alpha_search)*length(lambda_search),ncol=3)
# list_BIC[,1]=rep(alpha_search,each=length(lambda_search))
# list_BIC[,2]=rep(lambda_search,times=length(alpha_search))
# 
# for(aa in 1:nrow(list_BIC)){
#   list_BIC[aa,3]<-EM(y=y,k=max_k,alpha=list_BIC[aa,1],lambda=list_BIC[aa,2],size_factors<-size_factors)$BIC
#   print(list_BIC[aa,])
# }
# 
# max_index<-which(list_BIC[,3]==min(list_BIC[,3]))
# 
# max_alpha<-list_BIC[max_index,1]
# max_lambda<-list_BIC[max_index,2]
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
