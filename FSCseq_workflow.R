load("DESeq2")

source("FSCseq.R")

load("test_sim_set_K3_n150_g400.RData")     # create example_dat.RData
  # true.K = 3
  # log2 fold.change = 1
  # num.disc = 0.05
  # g = 2000
  # n = 150
  # distrib = "nb"
  # disp = "gene"
  # fixed_parms = T, fixed_coef = 6.5, fixed_phi = 0.15
  # X = matrix(batch,ncol=1), or X = NA (no covariates)

fit_DESeq_intercept=function(y,calc_vsd=F,calc_rld=F){
  n=ncol(y)
  coldata<-data.frame(matrix(rep(1,n),nrow=n))
  rownames(coldata)<-colnames(y)
  colnames(coldata)<-"int_only"
  dds<-DESeqDataSetFromMatrix(countData = y,
                              colData = coldata,
                              design = ~ 1)
  vsd=NA
  rld=NA
  size_factors=NA
  norm_y=NA
  
  dds=tryCatch({
    dds<-DESeq(dds)
  },error=function(err){
    tryCatch({
      print("All genes contain a 0. Using alternate size factor estimation")
      dds<-DESeq(dds, sfType = "iterate")
      return(dds)
    },error=function(err){
      print("Alternate SF estimation didn't converge. Trying zero-inflated model")
      dds<-DESeq(dds, sfType = "poscounts", useT=T, minmu=1e-6,minReplicatesForReplace = Inf)
      return(dds)
    })
  })
  
  if(!is.na(dds)){
    if(calc_vsd){
      vsd=varianceStabilizingTransformation(dds)
    }
    if(calc_rld){
      rld=rlogTransformation(dds)
    }
    size_factors<-sizeFactors(dds)
    norm_y<-counts(dds,normalized=TRUE)
  }
  
  return(list(dds=dds,vsd=vsd,rld=rld,size_factors=size_factors,norm_y=norm_y))
}

order_select = function(X=NA, y, size_factors, norm_y, K_search=c(1:7),             # since DESeq2 can take a long time for larger datasets, recommend it be run first
              distrib="nb", method="EM", disp="gene",
              true_clusters=NA, true_disc=NA,     # for diagnostics tracking ARI and nondisc/disc genes
              prefix="", dir="test"              # for diagnostics tracking file name/directory
              ){
  list_BIC=matrix(0,nrow=length(K_search),ncol=2)
  list_BIC[,1]=K
  list_res = list()
  
  for(aa in 1:nrow(list_BIC)){
    start=as.numeric(Sys.time())
    res<-FSCseq(X=X,y=y,k=list_BIC[aa,1],lambda=0,alpha=0,size_factors=size_factors,norm_y=norm_y,
                true_clusters=true_clusters,true_disc=true_disc,prefix=prefix,dir=dir,method=method,disp=disp)  # alpha = 0: all L1. No penalty here
    end=as.numeric(Sys.time())
    list_BIC[aa,2]<-res$BIC
    print(list_BIC[aa,])
    print(paste("Time:",end-start,"seconds"))
    list_res[[aa]]=res
  }
  
  names(list_res)=K
  max_k=list_BIC[which.min(list_BIC[,2]),1]
  unpen_BIC = min(list_BIC[,2])
  
  results = list(k = max_k,
                 BIC = list_BIC,
                 all_fits = list_res)
}

penalty_search = function(X=NA, y, size_factors, norm_y, 
                          distrib="nb", method="EM", disp="gene",
                          true_clusters=NA, true_disc=NA,
                          prefix="", dir="test",
                          lambda_search=c(0.1,seq(0.25,3,0.25)), alpha_search=seq(0,0.5,0.05),
                          OS_res = NA,                  # if there are OS results, then run that.
                          fit=NA, k=NA              # alternative to include custom k and unpenalized fit. must have FSCseq fit with coefs, phi, and clusters
                          ){
  
  if(exists("OS_res")){
    unpen_res = OS_res$all_fits[[which.min(list_BIC[,2])]]
    max_k = OS_res$k
  } else{
    unpen_res = fit
    max_k=k
  }
  
  list_BIC=matrix(0,nrow=length(lambda_search)*length(alpha_search),ncol=3) # matrix of BIC's: one for each combination of penalty params 
  
  list_BIC[,1]=rep(lambda_search,each=length(alpha_search))
  list_BIC[,2]=rep(alpha_search,times=length(lambda_search))
  
  list_pen_res = list()
  
  #search for optimal penalty parameters
  for(aa in 1:nrow(list_BIC)){
    pref = sprintf("grid%d",ii)
    start = as.numeric(Sys.time())
    res<-FSCseq(X=X,y=y,k=max_k,lambda=list_BIC[aa,1],alpha=list_BIC[aa,2],size_factors=size_factors,norm_y=norm_y,
                true_clusters=true_clusters,true_disc=true_disc,disp=disp,prefix=pref,dir=dir_name,method=method,
                init_cls=unpen_res$final_clusters,init_parms=T,init_coefs=unpen_res$coefs,init_phi=unpen_res$phi)
    end = as.numeric(Sys.time())
    list_BIC[aa,3]<-res$BIC
    print(list_BIC[aa,])
    print(paste("Time:",end-start,"seconds"))
    list_pen_res[[aa]]=res
  }
  
  #store optimal penalty parameters
  max_index<-which(list_BIC[,3]==min(list_BIC[,3]))
  max_lambda<-list_BIC[max_index,1]
  max_alpha<-list_BIC[max_index,2]
  
  print(paste("Dataset ", ii, "grid search results: ",list_BIC[max_index,]))
  
  if(length(max_index)>1){
    warning("more than one max index")
    max_index<-max_index[1]
    max_lambda<-list_BIC[max_index,1]
    max_alpha<-list_BIC[max_index,2]
  }
  
  pen_BIC = min(list_BIC[,3])
  
  if(pen_BIC<=unpen_res$BIC){
    final_res = list_pen_res[[max_index[1]]]
  } else{
    final_res = unpen_res
    print("Algorithm selected unpenalized model.")
  }
  
  results = list(fit = final_res,
                 k = max_k, cls = final_res$final_clusters,
                 lambda = max_lambda,
                 alpha = max_alpha,
                 all_fits = list_pen_res)
  
  return(results)
                 
}

# Calculate norm_y/SFs if not available:
fit=fit_DESeq_intercept(y,F,F)
norm_y=fit$norm_y
size_factors=fit$size_factors

# OS:
OS_res = order_select(X=NA, y=y, size_factors=size_factors, norm_y=norm_y, K_search=c(1:7),             # since DESeq2 can take a long time for larger datasets, recommend it be run first
                      distrib="nb", method="EM", disp="gene",
                      true_clusters=true_clusters, true_disc=true_disc,     # for diagnostics tracking ARI and nondisc/disc genes
                      prefix="", dir="test"              # for diagnostics tracking file name/directory)
                      )

PP_res = penalty_search(X=NA, y=y, size_factors=size_factors, norm_y=norm_y, 
                        distrib="nb", method="EM", disp="gene",
                        true_clusters=true_clusters, true_disc=true_disc,
                        prefix="", dir="test",
                        lambda_search=seq(0.05,1,0.05), alpha_search=seq(0,0.5,0.05),
                        OS_res = OS_res,                  # if there are OS results, then run that.
                        fit=NA, k=NA              # alternative to include custom k and unpenalized fit. must have FSCseq fit with coefs, phi, and clusters
                        )

final_res = PP_res$final_res
K = PP_res$k
cls = PP_res$cls
lambda = PP_res$lambda
alpha = PP_res$alpha

