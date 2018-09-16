
library(iClusterPlus)
library(DESeq2)
library(mclust)
library(parallel)


# DESeq analysis
normalizations <- function(dat){
  n=ncol(dat)
  g=nrow(dat)
  row_names<-paste("gene",seq(g))
  col_names<-paste("subj",seq(n))
  cts<-round(as.matrix(dat),digits=0)
  rownames(cts)<-row_names
  colnames(cts)<-col_names
  coldata<-data.frame(matrix(paste("cl"),nrow=n))
  rownames(coldata)<-colnames(cts)
  colnames(coldata)<-"cluster"
  dds<-DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ 1)
  dds<-DESeq(dds)
  size_factors<-sizeFactors(dds)
  norm_y<-counts(dds,normalized=TRUE)
  idx <- rowSums( counts(dds, normalized=TRUE) >= 20) >= 5
  
  rownames(norm_y) = rownames(dat)
  colnames(norm_y) = colnames(dat)
  
  res <- list(size_factors=size_factors,
              norm_y=norm_y,
              idx=idx)
  return(res)
}




compare <- function(y){
  
  true_cls = as.numeric(NSCLC_anno$Adeno.Squamous)
  true_k = 4
  
  # Number of cores: 2 for laptop, 12 for Killdevil
  #no_cores <- 2   # for parallel computing
  no_cores <- 10
  
  # list of K to search over
  K_search = c(2:7)

  
  # iClusterPlus #
  i_tune_start_time <- Sys.time()
  
  iClust_OS <- list()
  for(k in (K_search-1)){
    cv.fit = tune.iClusterPlus(cpus=no_cores,dt1=t(y),type="poisson",K=k,alpha=1,n.lambda=25,scale.lambda=1,maxiter=20)
    iClust_OS[[k]] = cv.fit
  }
  BIC_mat = getBIC(iClust_OS)
  
  BIC_OS = apply(BIC_mat,2,min)
  lambda_OS = iClust_OS[[1]]$lambda[apply(BIC_mat,2,which.min),1]
  
  i_max_k = which.min(BIC_OS)
  max_lambda = lambda_OS[i_max_k]
  
  i_tune_end_time <- Sys.time()
  i_tune_time <- as.numeric(i_tune_end_time)-as.numeric(i_tune_start_time)
  
  i_start <- Sys.time()
  iClust_fit <- iClusterPlus(dt1=t(y),type="poisson",lambda=max_lambda,alpha=1,K=i_max_k,maxiter=10)   # K=1 is true
  i_end <- Sys.time()
  i_time <- as.numeric(i_end)-as.numeric(i_start)
  
  icoefs <- cbind(iClust_fit$alpha[[1]],iClust_fit$alpha[[1]]+iClust_fit$beta[[1]])
  
  
  features = rownames(y)
  
  rowsum=apply(abs(iClust_fit$beta[[1]]),1, sum)
  upper=quantile(rowsum,prob=0.75)
  isigfeatures=(features)[which(rowsum>upper)]
  
  print(paste("Order:",(i_max_k+1)))
  print(paste("Clusters:",iClust_fit$clusters))
  
  
  
  
  # my method #
  source("NB Pan EM par.R")
  dir_name = "BRCA"
  
  EM_tune_start <- Sys.time()
  
  list_BIC=matrix(0,nrow=length(K_search),ncol=2)
  list_BIC[,1]=K_search
  
  
  # for(aa in 1:nrow(list_BIC)){
  #   list_BIC[aa,2]<-EM(y=y,k=list_BIC[aa,1],lambda1=0,lambda2=0,tau=0,size_factors=size_factors,norm_y=norm_y)$BIC   # no penalty Pan
  #   print(list_BIC[aa,])
  # }
  
  order.select <- function(c){
    X<-EM(y=y,k=c,lambda1=0,lambda2=0,tau=0,size_factors=size_factors,norm_y=norm_y,true_clusters=true_cls,prefix="order",dir=dir_name)
    return(X$BIC)
  }
  
  cl<-makeCluster(no_cores)
  clusterExport(cl=cl,varlist=c(ls(envir=environment()),"y","EM","EM_run","logsumexpc","soft_thresholding"),envir=environment())
  clusterExport(cl=cl,varlist=c(ls(envir=globalenv()),"norm_y","size_factors"),envir=globalenv())
  clusterEvalQ(cl,{
    library(stats)
    library(MASS)
    library(permute)
    library(Rcpp)
    library(RcppArmadillo)
    library(mclust)
    library(pryr)
    sourceCpp("M_step.cpp")
  })
  
  list_BIC[,2] <- parSapply(cl,list_BIC[,1],order.select)
  stopCluster(cl)
  
  max_k=list_BIC[which(list_BIC[,2]==min(list_BIC[,2])),1]
  
  fit = EM(y=y,k=c,lambda1=0,lambda2=0,tau=0,size_factors=size_factors,norm_y=norm_y,true_clusters=true_cls,dir=dir_name)
  compare_finalcls = fit$final_clusters
  compare_wts = fit$wts
  
  sink(file=sprintf("Diagnostics/%s/final_order.txt",dir_name))
  
  max_k=list_BIC[which.min(list_BIC[,2]),1]
  print(paste("True order:",true_k))
  print(paste("Optimal order selected:",max_k))
  print("RUN WITH CORRECT ORDER:")
  print(paste("ARI:",adjustedRandIndex(compare_finalcls,true_cls)))
  print(paste("PPs:"))
  print(compare_wts)
  
  sink()
  
  pdf(file=sprintf("Diagnostics/%s/final_order.pdf",dir_name))
  for(c in 1:true_k){
    cl_ids = true_cls==c
    boxplot(compare_wts[c,cl_ids],main=sprintf("Boxplot of PP for cl%d",c))
  }
  dev.off()
  
  
  
  
  lambda1_search=seq(from=0.1,to=2.1,by=0.2)
  lambda2_search=c(0.5,1,1.5,2,2.5,3,3.5,4)
  tau_search=seq(from=0.1,to=0.9,by=0.2)
  
  list_BIC=matrix(0,nrow=length(lambda1_search)*length(lambda2_search)*length(tau_search),ncol=4) # matrix of BIC's: one for each combination of penalty params 
  list_BIC[,1]=rep(lambda1_search,each=length(lambda2_search)*length(tau_search))
  list_BIC[,2]=rep(rep(lambda2_search,each=length(tau_search)),times=length(lambda1_search))
  list_BIC[,3]=rep(tau_search,times=length(lambda1_search)*length(lambda2_search))
  
  #search for optimal penalty parameters
  # for(aa in 1:nrow(list_BIC)){
  #   X<-EM(y=y,k=max_k,tau=list_BIC[aa,3],lambda1=list_BIC[aa,1],lambda2=list_BIC[aa,2],size_factors=size_factors,norm_y=norm_y,true_clusters=NA)
  #   list_BIC[aa,4]<-X$BIC
  #   print(list_BIC[aa,])
  #   print(paste("Time:",X$time_elap,"seconds"))
  # }
  
  grid.search <- function(row){
    X<-EM(y=y,k=max_k,lambda1=list_BIC[row,3],lambda2=list_BIC[row,1],tau=list_BIC[row,2],size_factors=size_factors,norm_y=norm_y,true_clusters=true_cls,prefix="tuning",dir=dir_name)
    return(X$BIC)
  }
  cl<-makeCluster(no_cores)
  clusterExport(cl=cl,varlist=c(ls(envir=environment()),"y","max_k","EM","EM_run","logsumexpc","soft_thresholding"),envir=environment())
  clusterExport(cl=cl,varlist=c(ls(envir=globalenv()),"norm_y","size_factors"),envir=globalenv())
  clusterEvalQ(cl,{
    library(stats)
    library(MASS)
    library(permute)
    library(Rcpp)
    library(RcppArmadillo)
    library(mclust)
    library(pryr)
    sourceCpp("M_step.cpp")
  })
  
  list_BIC[,4] <- parSapply(cl,c(1:nrow(list_BIC)),grid.search)
  stopCluster(cl)
  
  max_index<-which(list_BIC[,4]==min(list_BIC[,4]))
  max_tau<-list_BIC[max_index,3]
  max_lambda1<-list_BIC[max_index,1]
  max_lambda2<-list_BIC[max_index,2]
  
  print(list_BIC[max_index,])
  
  if(length(max_index)>1){
    warning("more than one max index")
    max_index<-max_index[1]
    max_tau<-list_BIC[max_index,3]
    max_lambda1<-list_BIC[max_index,1]
    max_lambda2<-list_BIC[max_index,2]
  }
  
  EM_tune_end <- Sys.time()
  
  EM_tune_time <- as.numeric(EM_tune_end)-as.numeric(EM_tune_start)
  
  EM_start <- Sys.time()
  X<-EM(y=y,k=max_k,tau=max_tau,lambda1=max_lambda1,lambda2=max_lambda2,size_factors=size_factors,norm_y=norm_y,true_clusters=true_cls,prefix="final",dir=dir_name)
  EM_end <- Sys.time()
  EM_time <- as.numeric(EM_end)-as.numeric(EM_start)
  
  EM_sigfeatures<-features[X$nondiscriminatory==FALSE]
  
  list(i_maxK=(i_max_k+1),
       EM_maxK=max_k,
       i_sigfeatures=isigfeatures,
       EM_sigfeatures=EM_sigfeatures,
       i_clusters=iClust_fit$clusters,
       EM_clusters=X$final_clusters,
       i_tune_time=i_tune_time,
       i_time=i_time,
       EM_tune_time=EM_tune_time,
       EM_time=EM_time,
       i_X=iClust_fit,
       EM_X=X)
}

#setwd("/Users/deelim/Documents/Research")
setwd("/netscr/deelim/out")

# ################################ NSCLC (Cell Line) Data #######################################
# 
# # DATA PREP
# 
# NSCLC_anno <- read.table("NSCLC_anno.txt",sep="\t",header=TRUE)
# NSCLC_dat <- read.table("NSCLC_rsem.genes.exp.count.unnormalized.txt",sep="\t",header=TRUE)
# rownames(NSCLC_dat)<-toupper(NSCLC_dat[,1])
# NSCLC_subtypes <- NSCLC_anno$Adeno.Squamous
# NSCLC_dat<-round(NSCLC_dat[,-1],digits=0)
# 
# 
# # pre-filtering DESEQ significant genes
# 
# ## DESeq2 to find top significant genes in clustering
# # colnames(NSCLC_dat)<-toupper(colnames(NSCLC_dat))
# # coldata<-NSCLC_anno[,-1]
# # rownames(coldata)<-toupper(NSCLC_anno[,1])
# # coldata<-coldata[,c("Adeno.Squamous","Tumor.location")]
# # dds<-DESeqDataSetFromMatrix(countData = NSCLC_dat,
# #                             colData = coldata,
# #                             design = ~ Adeno.Squamous)
# # DESeq_dds<-DESeq(dds)
# 
# 
# # RUN
# fit <- normalizations(NSCLC_dat)    # Normalizations done on full dataset
# size_factors <- fit$size_factors
# norm_y <- fit$norm_y
# 
# filt_id = fit$idx & rowMedians(norm_y)>=15           # pre-filtering rowmedians >=100
# #rm(fit)
# 
# y1 <- NSCLC_dat[filt_id,]
# norm_y = norm_y[filt_id,]
# rownames(y1) = rownames(NSCLC_dat)[filt_id]
# rownames(norm_y) = rownames(y1)
# colnames(y1) = colnames(NSCLC_dat)
# colnames(norm_y) = colnames(y1)
# 
# X1<-compare(y1)
# 






################################ BRCA (Breast Cancer) Data #######################################

# library(SummarizedExperiment)
# library(genefilter)
#
# # DATA PREP
#
# load(file="TCGA_BRCA_exp.rda")
# BRCA_cts <- round(assay(data),digits=0)             # default gives raw count
# BRCA_anno <- colData(data)@listData
# rm(data)
#
# BRCA_subtypes <- BRCA_anno$subtype_PAM50.mRNA
#
# idy <- (!is.na(BRCA_subtypes) & BRCA_subtypes != "Normal-like")
# BRCA_cts <- BRCA_cts[!duplicated(BRCA_cts[,1:ncol(BRCA_cts)]),idy]
# BRCA_subtypes = BRCA_subtypes[idy]




# RUN


X2 <- compare(y2)


save(X2,file="BRCA_compare.out")



# 
# load("X1.out")
# load("X2.out")
# 
# rowsum1=apply(abs(X1$i_X$beta[[1]]),1,sum)
# cutoff1=(X1$EM_X$lambda2)
# features1=rownames(y1)
# isigfeatures1=features1[which(rowsum1>cutoff1)]
# 
# rowsum2=apply(abs(X2$i_X$beta[[1]]),1, sum)
# cutoff2=X2$EM_X$lambda2
# features2=rownames(y2)
# isigfeatures2=(features2)[which(rowsum2>cutoff2)]
# 
# write_csv = function(X,y){
#   icoefs = cbind(X$i_X$alpha[[1]], X$i_X$alpha[[1]] + X$i_X$beta[[1]])
#   EMcoefs = X$EM_X$coefs
#   
#   i_FC = (apply(icoefs,1,max) - apply(icoefs,1,min))/(X$i_maxK - 1)      # avg FC = (Max - min) / (#clusters - 1)
#   EM_FC = (apply(EMcoefs,1,max) - apply(EMcoefs,1,min))/(X$EM_maxK - 1)
#   
#   lambda2 = X$EM_X$lambda2
#   features=rownames(y)
#   
#   isig_id=which(i_FC>lambda2)
#   EMsig_id1=(X$EM_X$nondiscriminatory==FALSE) # EM theta = 0 for all K as criteria
#   EMsig_id2=which(EM_FC>lambda2)       # EM FC > lambda2 as criteria for discriminatory
#                                              # These criteria are identical for K=2, but not in general
#   isig_feat=features[isig_id]
#   EMsig_feat1=features[EMsig_id1]
#   EMsig_feat2=features[EMsig_id2]
#   
#   
#   # Order by FC
#   i_order = order(-i_FC)
#   EM_order = order(-EM_FC)
#   
#   tab <- cbind(features[i_order],i_FC[i_order],(features[i_order] %in% isig_feat)^2,
#                features[EM_order],EM_FC[EM_order],(features[EM_order] %in% EMsig_feat1)^2,(features[EM_order] %in% EMsig_feat2)^2)
#   colnames(tab) <- c("i_Features","i_FC","i_sig","EM_Features","EM_FC","EM_sig_theta","EM_sig_FC")
# 
#   return(tab)
# }
# 
# tab1<-write_csv(X1,y1)
# tab2<-write_csv(X2,y2)
# 
# write.csv(tab1,"NSCLC_comparison.csv")
# write.csv(tab2,"BRCA_comparison.csv")


################################################
#d<-as.dist(1-cor(norm_y, method="spearman"))  ##Spearman correlation distance w/ log transform##
#model<-hclust(d,method="complete")       # hierarchical clustering
#heatmap(as.matrix(norm_y),Rowv = as.dendrogram(model),Colv=as.dendrogram(model))
###############################################