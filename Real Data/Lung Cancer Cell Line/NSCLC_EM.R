#setwd("C:/Users/David/Desktop/Research/GitHub/EM/Real Data/Lung Cancer Cell Line")
setwd("/netscr/deelim")
#source("C:/Users/David/Desktop/Research/GitHub/EM/non-simulation EM.R")
source("non-simulation EM.R")
library("parallel")

no_cores<-detectCores()-1
cl<-makeCluster(no_cores,outfile="NSCLC_EM_cluster_output.txt")

library("stats")
library("data.table")
library("DESeq2")

anno<-read.table("NSCLC_anno.txt",sep="\t",header=TRUE)
dat<-read.table("NSCLC_rsem.genes.exp.count.unnormalized.txt",sep="\t",header=TRUE)
rownames(cts)<-toupper(dat[,1])
dat<-round(dat[,-1],digits=0)

# DESeq2 to find size factors
cts<-as.matrix(dat)
colnames(cts)<-toupper(colnames(cts))
coldata<-anno[,-1]
rownames(coldata)<-toupper(anno[,1])
coldata<-coldata[,c("Adeno.Squamous","Tumor.location")]
all(rownames(coldata) %in% colnames(cts))
all(rownames(coldata) == colnames(cts))
dds<-DESeqDataSetFromMatrix(countData = cts,
                            colData = coldata,
                            design = ~ 1)
DESeq_dds<-DESeq(dds)
size_factors<-sizeFactors(DESeq_dds)


# initial clean-up of data and pre-filtering to include only genes with >=100 count
y<-dat[,-1]
rownames(y)<-dat[,1]
y<-y[(rowSums(y)>=100),]


# filtering genes to have top 500 MAD (median absolute deviation): optional
#med_abs_dev<-rep(0,times=nrow(y))
#for(j in 1:nrow(y)){
#  med_abs_dev[j]<-mad(as.numeric(y[j,]),constant=1)
#}
#y<-cbind(rownames(y),y,med_abs_dev)
#subs_y<-as.data.table(y)[order(-med_abs_dev),head(.SD,500)]
#genes_y<-subs_y[,1]
#subs_y<-subs_y[,-1]
#subs_y<-as.data.frame(subs_y[,-ncol(subs_y)])



# k=2        # known

# grid search for tuning params lambda1 and lambda2 and K

lambda1_search=c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1)
lambda2_search=c(0.01, 0.05, 0.1, 0.2, 0.5, 0.8, 1)
K_search=c(2:6)

list_BIC=matrix(0,nrow=length(lambda1_search)*length(lambda2_search)*length(K_search),ncol=4) #matrix of BIC's: lambda1 and lambda2 and K, 49*5 combinations

list_BIC[,1]=rep(lambda1_search,each=length(lambda2_search)*length(K_search))
list_BIC[,2]=rep(lambda2_search,times=length(lambda1_search)*length(K_search))
list_BIC[,3]=rep(rep(K_search,each=length(lambda1_search)),times=length(lambda2_search))

clusterExport(cl,c("y","size_factors","list_BIC","EM","logsumexpc","soft_thresholding"))

extract_BIC<-function(row){
  print(paste("lambda1 =",list_BIC[row,1],"and lambda2 =",list_BIC[row,2],"and K =",list_BIC[row,3]))
  X<-EM(y=y,k=list_BIC[row,3],lambda1=list_BIC[row,1],lambda2=list_BIC[row,2],size_factors=size_factors)
  print(X$pi)
  print(X$BIC)
  print(X$nondiscriminatory)
  return(X$BIC)
}

clusterExport(cl,"extract_BIC")


# actual grid search run
list_BIC[,4]<-parSapply(cl,1:nrow(list_BIC),extract_BIC)


# storing optimal BIC index & optimal parameters
max_index<-which(list_BIC[,4]==min(list_BIC[,4]))

max_lambda1<-list_BIC[max_index,1]
max_lambda2<-list_BIC[max_index,2]
max_k<-list_BIC[max_index,3]


# actual run:
X<-EM(y=y,k=max_k,lambda1=max_lambda1,lambda2=max_lambda2,size_factors=size_factors)


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