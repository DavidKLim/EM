#setwd("C:/Users/David/Desktop/Research/GitHub/EM/Real Data/Lung Cancer Cell Line")
setwd("/netscr/deelim")
#source("C:/Users/David/Desktop/Research/GitHub/EM/non-simulation EM.R")
source("non-simulation EM.R")
library("parallel")

no_cores<-detectCores()-1
cl<-makeCluster(no_cores)

library("stats")
library("data.table")

anno<-read.table("NSCLC_anno.txt",sep="\t",header=TRUE)
dat<-read.table("NSCLC_rsem.genes.exp.count.unnormalized.txt",sep="\t",header=TRUE)
y<-dat[,-1]
rownames(y)<-dat[,1]
y<-y[(rowSums(y)>=100),]
y<-round(y,digits=0)

med_abs_dev<-rep(0,times=nrow(y))
for(j in 1:nrow(y)){
  med_abs_dev[j]<-mad(as.numeric(y[j,]),constant=1)
}

y<-cbind(rownames(y),y,med_abs_dev)

subs_y<-as.data.table(y)[order(-med_abs_dev),head(.SD,500)]

genes_y<-subs_y[,1]
subs_y<-subs_y[,-1]

subs_y<-as.data.frame(subs_y)


# k=2        # known

# search for optimal k using BIC
clusterExport(cl,c("subs_y","EM","logsumexpc","soft_thresholding"))

#for(k in 2:6){
#  k_BIC[k]<-EM(y=y,k=k,lambda1=0.5,lambda2=0.5)$BIC
#}

#k_BIC<-rep(0,times=5)
#k_BIC<-parSapply(cl,2:6,function(k) EM(y=subs_y,k=k,lambda1=0.5,lambda2=0.5)$BIC)
#max_k<-which.max(k_BIC)+1
#print(max_k)


# grid search for tuning params lambda1 and lambda2

lambda1_search=c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1)
lambda2_search=c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1)
K_search=c(2:6)

list_BIC=matrix(0,nrow=length(lambda1_search)*length(lambda2_search)*length(K_search),ncol=4) #matrix of BIC's: lambda1 and 2, 49 combinations

list_BIC[,1]=rep(lambda1_search,each=length(lambda2_search)*length(K_search))
list_BIC[,2]=rep(lambda2_search,times=length(lambda1_search)*length(K_search))
list_BIC[,3]=rep(rep(K_search,each=length(lambda1_search)),times=length(lambda2_search))

extract_BIC<-function(row){
  print(paste("lambda1 =",list_BIC[row,1],"and lambda2 =",list_BIC[row,2],"and K =",list_BIC[row,3]))
  X<-EM(y=subs_y,k=list_BIC[row,3],lambda1=list_BIC[row,1],lambda2=list_BIC[row,2])
  print(X$pi)
  print(X$BIC)
  print(X$nondiscriminatory)
  return(X$BIC)
}


clusterExport(cl,c("list_BIC","extract_BIC"))


#for(aa in 1:nrow(list_BIC)){
#  set.seed(aa)
#  list_BIC[aa,3]<-EM(y=y,k=max_k,lambda1=list_BIC[aa,1],lambda2=list_BIC[aa,2])$BIC
#  print(list_BIC[aa,])
#}
list_BIC[,4]<-parSapply(cl,1:nrow(list_BIC),extract_BIC)

max_index<-which(list_BIC[,4]==min(list_BIC[,4]))

max_lambda1<-list_BIC[max_index,1]
max_lambda2<-list_BIC[max_index,2]
max_k<-list_BIC[max_index,3]


# actual run:
X<-EM(y=subs_y,k=max_k,lambda1=max_lambda1,lambda2=max_lambda2)


stopCluster(cl)

# summarize output #
sink("NSCLC_EM.txt",append=FALSE)
print("Results")
print(paste("pi =",X$pi))
print(paste("mean % of nondiscriminatory genes =",X$nondiscriminatory))
print(paste("final (lambda1,lambda2) =",max_lambda1,max_lambda2))
print(paste("final clusters:",X$final_clusters))
sink()

sink("NSCLC_EM_coefs.txt",append=FALSE)
print("coefs")
print(X$coefs)
sink()