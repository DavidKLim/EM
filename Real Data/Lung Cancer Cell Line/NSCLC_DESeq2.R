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