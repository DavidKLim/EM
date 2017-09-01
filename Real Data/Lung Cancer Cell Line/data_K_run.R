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
y<-y+1


source("Pan EM.R")

K_search=c(2:7)
list_BIC=matrix(0,nrow=length(K_search),ncol=2)
list_BIC[,1]=K_search

extract_BIC<-function(row){
  X<-EM(y=y,k=list_BIC[row,1],lambda1=0,lambda2=0,tau=9999,size_factors=size_factors)
  print(paste("K =",list_BIC[row,1],"pi=",X$pi,"BIC=",X$BIC,"nondisc=",mean(X$nondiscriminatory)))
  return(X$BIC)
}
clusterExport(cl,c("y","size_factors","norm_y","list_BIC","EM","logsumexpc","soft_thresholding"))
clusterExport(cl,"extract_BIC")

list_BIC[,2]<-parSapply(cl,1:nrow(list_BIC),extract_BIC)



stopCluster(cl)

