setwd("C:/Users/David/Desktop/Research/GitHub/EM/Real Data/Lung Cancer Cell Line")
setwd("/netscr/deelim")
library("parallel")

no_cores<-detectCores()-1
cl<-makeCluster(no_cores,outfile="NSCLC_EM_cluster_output.txt")

library("stats")
library("data.table")
library("DESeq2")

anno<-read.table("NSCLC_anno.txt",sep="\t",header=TRUE)
dat<-read.table("NSCLC_rsem.genes.exp.count.unnormalized.txt",sep="\t",header=TRUE)
row_names<-toupper(dat[,1])
dat<-round(dat[,-1],digits=0)

# DESeq2 to find size factors
cts<-as.matrix(dat)
rownames(cts)<-row_names

colnames(cts)<-toupper(colnames(cts))
coldata<-anno[,-1]
rownames(coldata)<-toupper(anno[,1])
coldata<-coldata[,c("Adeno.Squamous","Tumor.location")]
#all(rownames(coldata) %in% colnames(cts))         # check that headers are correct
#all(rownames(coldata) == colnames(cts))
dds<-DESeqDataSetFromMatrix(countData = cts,
                            colData = coldata,
                            design = ~ 1)
DESeq_dds<-DESeq(dds)
size_factors<-sizeFactors(DESeq_dds)

norm_y<-counts(DESeq_dds,normalized=TRUE)



# initial clean-up of data and pre-filtering to include only genes with >=100 count
dat<-read.table("NSCLC_rsem.genes.exp.count.unnormalized.txt",sep="\t",header=TRUE)
y<-round(dat[,-1],digits=0)
rownames(y)<-dat[,1]
y<-y[(rowSums(y)>=100),]


# filtering genes to have top 500 MAD (median absolute deviation): optional
med_abs_dev<-rep(0,times=nrow(y))
for(j in 1:nrow(y)){
  med_abs_dev[j]<-mad(as.numeric(y[j,]),constant=1)
}
y<-cbind(rownames(y),y,med_abs_dev)
subs_y<-as.data.table(y)[order(-med_abs_dev),head(.SD,1000)]
genes_y<-subs_y[,1]
subs_y<-subs_y[,-1]
subs_y<-as.data.frame(subs_y[,-24])



# k=2        # known

# grid search for tuning params lambda1 and lambda2 and K
# Wei Pan
#source("C:/Users/David/Desktop/Research/GitHub/EM/Pan EM.R")
source("Pan EM.R")
lambda1_search=c(0.01,0.05,0.1,0.2,1,1.5,2)
lambda2_search=c(0.01,0.05,0.1,0.2,1,1.5,2,5,10,100,1000)
tau_search=seq(from=0.1,to=1,by=0.1)
K_search=c(2:6)

list_BIC=matrix(0,nrow=length(lambda1_search)*length(lambda2_search)*length(K_search)*length(tau_search),ncol=5) #matrix of BIC's: lambda1 and lambda2 and K, 49*5 combinations

list_BIC[,1]=rep(lambda1_search,each=length(lambda2_search)*length(K_search)*length(tau_search))
list_BIC[,2]=rep(rep(lambda2_search,each=length(K_search)*length(tau_search)),times=length(lambda1_search))
list_BIC[,3]=rep(rep(K_search,each=length(tau_search)),times=length(lambda1_search)*length(lambda2_search))
list_BIC[,4]=rep(tau_search,times=length(lambda1_search)*length(lambda2_search)*length(K_search))

extract_BIC<-function(row){
  X<-EM(y=subs_y,k=list_BIC[row,3],lambda1=list_BIC[row,1],lambda2=list_BIC[row,2],tau=list_BIC[row,4],size_factors=size_factors)
  print(paste("lambda1 =",list_BIC[row,1],"and lambda2 =",list_BIC[row,2],"and K =",list_BIC[row,3],"and tau =",list_BIC[row,4],"pi=",X$pi,"BIC=",X$BIC,"nondisc=",mean(X$nondiscriminatory)))
  return(X$BIC)
}
clusterExport(cl,c("subs_y","y","size_factors","norm_y","list_BIC","EM","logsumexpc","soft_thresholding"))
clusterExport(cl,"extract_BIC")
# actual grid search run
list_BIC[,5]<-parSapply(cl,1:nrow(list_BIC),extract_BIC) # Wei Pan

# storing optimal BIC index & optimal parameters
max_index<-which(list_BIC[,ncol(list_BIC)]==min(list_BIC[,ncol(list_BIC)]))
if(length(max_index)>1){
  warning("More than one max index")
  max_index<-max_index[1]
}

max_lambda1<-list_BIC[max_index,1]
max_lambda2<-list_BIC[max_index,2]
max_k<-list_BIC[max_index,3]
max_tau<-list_BIC[max_index,4]             # Wei Pan

print(paste("lambda1, lambda2, k, tau =", max_lambda1, max_lambda2, max_k, max_tau))    # Wei Pan

# actual run:
X_pan<-EM(y=subs_y,k=max_k,lambda1=max_lambda1,lambda2=max_lambda2,tau=max_tau,size_factors=size_factors) # Wei Pan


stopCluster(cl)

# summarize output #
sink("NSCLC_EM_Pan.txt",append=FALSE)
print("Results")
print(paste("pi =",X_pan$pi))
print(paste("mean % of nondiscriminatory genes =",X_pan$nondiscriminatory))
print(paste("final (lambda1,lambda2) =",max_lambda1,max_lambda2))   # Wei Pan
print(paste("final clusters:",X_pan$final_clusters))
sink()

sink("NSCLC_EM_coefs_Pan.txt",append=FALSE)
print("coefs")
print(X_pan$coefs)
sink()
















# group lasso
#source("C:/Users/David/Desktop/Research/GitHub/EM/non-simulation EM group lasso.R")
source("non-simulation EM group lasso.R")
K_search=c(2:6)
alpha_search=seq(from=0,to=1,by=0.1)
lambda_search=seq(from=0,to=1,by=.2)
list_BIC=matrix(0,nrow=length(K_search)*length(alpha_search)*length(lambda_search),ncol=4)
list_BIC[,1]=rep(K_search,each=length(alpha_search)*length(lambda_search))
list_BIC[,2]=rep(rep(alpha_search,each=length(lambda_search)),times=length(K_search))
list_BIC[,3]=rep(lambda_search,times=length(K_search)*length(alpha_search))

extract_BIC<-function(row){
  X<-EM(y=subs_y,k=list_BIC[row,1],alpha=list_BIC[row,2],lambda=list_BIC[row,3],size_factors=size_factors)
  print(paste("k=",list_BIC[row,1],"alpha=",list_BIC[row,2],"lambda=",lambda=list_BIC[row,3],"BIC=",X$BIC))
  return(X$BIC)
}



clusterExport(cl,c("subs_y","y","size_factors","norm_y","list_BIC","EM","logsumexpc","soft_thresholding"))
clusterExport(cl,"extract_BIC")


# actual grid search run
list_BIC[,4]<-parSapply(cl,1:nrow(list_BIC),extract_BIC) # Group lasso


# storing optimal BIC index & optimal parameters
max_index<-which(list_BIC[,ncol(list_BIC)]==min(list_BIC[,ncol(list_BIC)]))
if(length(max_index)>1){
  warning("More than one max index")
  max_index<-max_index[1]
}

max_k<-list_BIC[max_index,1]
max_alpha<-list_BIC[max_index,2]             
max_lambda<-list_BIC[max_index,3]            # Group lasso

print(paste("k=",max_k,"alpha=",max_alpha,"lambda=",max_lambda))          # Group lasso

# actual run:
X_glmnet<-EM(y=subs_y,k=max_k,alpha=max_alpha,lambda=max_lambda,size_factors=size_factors)     # Group lasso

stopCluster(cl)

# summarize output #
sink("NSCLC_EM_glmnet.txt",append=FALSE)
print("Results")
print(paste("pi =",X_glmnet$pi))
print(paste("mean % of nondiscriminatory genes =",X_glmnet$nondiscriminatory))
print(paste("final alpha,lambda =",max_alpha,max_lambda))                              # Group lasso
print(paste("final clusters:",X_glmnet$final_clusters))
sink()

sink("NSCLC_EM_coefs_glmnet.txt",append=FALSE)
print("coefs")
print(X_glmnet$coefs)
sink()