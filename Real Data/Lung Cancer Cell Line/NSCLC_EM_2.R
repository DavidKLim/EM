#setwd("C:/Users/David/Desktop/Research/GitHub/EM/Real Data/Lung Cancer Cell Line")
setwd("/netscr/deelim")
#source("C:/Users/David/Desktop/Research/GitHub/EM/non-simulation EM.R")
source("non-simulation EM.R")

library("parallel")
no_cores<-detectCores()-1
cl<-makeCluster(no_cores)



anno<-read.table("NSCLC_anno.txt",sep="\t",header=TRUE)
dat<-read.table("NSCLC_rsem.genes.exp.count.unnormalized.txt",sep="\t",header=TRUE)
y<-dat[,-1]
y<-y[(rowSums(y)>=100),]
y<-round(y,digits=0)
# k=2        # known

# search for optimal k using BIC
clusterExport(cl,c("y","EM","logsumexpc","soft_thresholding"))

k_BIC<-rep(0,times=5)
#for(k in 2:6){
#  k_BIC[k]<-EM(y=y,k=k,lambda1=0.5,lambda2=0.5)$BIC
#}
k_BIC<-parSapply(cl,2:6,function(k) EM(y=y,k=k,lambda1=0.5,lambda2=0.5)$BIC)
max_k<-which.max(k_BIC)+1
print(max_k)


# grid search for tuning params lambda1 and lambda2

lambda1_search=c(0.01,0.05, 0.1, 0.2, 0.5, 0.8, 1)
lambda2_search=c(0.01,0.05, 0.1, 0.2, 0.5, 0.8, 1)

list_BIC=matrix(0,nrow=length(lambda1_search)*length(lambda2_search),ncol=3) #matrix of BIC's: lambda1 and 2, 49 combinations

list_BIC[,1]=rep(lambda1_search,each=length(lambda2_search))
list_BIC[,2]=rep(lambda2_search,times=length(lambda1_search))

extract_BIC<-function(row){
  return(EM(y=y,k=max_k,lambda1=list_BIC[row,1],lambda2=list_BIC[row,2])$BIC)
}

clusterExport(cl,c("list_BIC","max_k","extract_BIC"))

#for(aa in 1:nrow(list_BIC)){
#  set.seed(aa)
#  list_BIC[aa,3]<-EM(y=y,k=max_k,lambda1=list_BIC[aa,1],lambda2=list_BIC[aa,2])$BIC
#  print(list_BIC[aa,])
#}
list_BIC[,3]<-parSapply(cl,1:nrow(list_BIC),extract_BIC)

max_index<-which(list_BIC[,3]==min(list_BIC[,3]))

max_lambda1<-list_BIC[max_index,1]
max_lambda2<-list_BIC[max_index,2]


# actual run:
clusterExport(cl,c("max_lambda1","max_lambda2"))
X<-EM(y=y,k=max_k,lambda1=max_lambda1,lambda2=max_lambda2)





stopCluster(cl)

# summarize output #
sink("NSCLC_EM.txt",append=FALSE)
print("Results")
print(paste("pi =",X$pi))
print(paste("mean % of nondiscriminatory genes =",X$nondiscriminatory))
print(paste("final (lambda1,lambda2) =",max_lambda1,max_lambda2))
sink()

sink("NSCLC_EM_coefs.txt",append=FALSE)
print("coefs")
print(X$coefs)
sink()