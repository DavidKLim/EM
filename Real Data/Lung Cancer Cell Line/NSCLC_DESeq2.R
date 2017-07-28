#setwd("C:/Users/David/Desktop/Research/GitHub/EM/Real Data/Lung Cancer Cell Line")
setwd("/netscr/deelim")

library("DESeq2")


anno<-read.table("NSCLC_anno.txt",sep="\t",header=TRUE)
dat<-read.table("NSCLC_rsem.genes.exp.count.unnormalized.txt",sep="\t",header=TRUE)
rownames(cts)<-toupper(dat[,1])
dat<-round(dat[,-1],digits=0)
# k=2        # known

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
