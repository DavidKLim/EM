#setwd("C:/Users/David/Desktop/Research/EM")

# NSCLC
#init_y<-read.table("init_y.txt")
#init_size_factors<-as.numeric(read.table("init_size_factors.txt")[,1])
#init_norm_y<-read.table("init_norm_y.txt")
#init_y<-cbind(init_y,init_y,init_y,init_y,init_y,init_y,init_y,init_y,init_y,init_y,init_y,init_y)

# BRCA
load('sim_y.RData')
init_y = cbind(y,y,y)
colnames(init_y) = 1:ncol(init_y)

library("stats")
library("data.table")
library("DESeq2")
library("mclust")
library("parallel")
library("pheatmap")
library("iClusterPlus")
library("cluster")
library("NbClust")
library("NB.MClust")

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

# Function to simulate data
simulate_data=function(n,k,g,init_pi,b,size_factors,distrib,phi=matrix(0,nrow=g,ncol=k),
                       batch_coef=1,batch=rep(0,n),batch_g=NA){      # batch of 0: no batch effect, batch of 1: yes effect
  batch_eff = batch_coef*batch
  y<-matrix(rep(0,times=g*n),nrow=g)
  z = rmultinom(n,1,init_pi)
  if(all(is.na(batch_g))){
    batch_g = sample(1:g, 0.5*g)   # Apply batch on 50% of genes
  }
  if(ncol(b)!=k){
    warning("Wrong order selected. Simulating based on correct order")
    k=ncol(b)
  }
  cl = rep(0,n)
  for(c in 1:k){
    cl[z[c,]==1] = c
  }
  if(distrib=="poisson"){
    for(j in 1:g){
      for(i in 1:n){
        if(j %in% batch_g){
          y[j,i] = rpois( 1, lambda = size_factors[i]*2^(b[j,cl[i]] + batch_eff[i]))
        } else{
          y[j,i] = rpois( 1, lambda = size_factors[i]*2^(b[j,cl[i]]))
        }
      }
    }
  } else if(distrib=="nb"){
    for(j in 1:g){
      for(i in 1:n){
        if(j %in% batch_g){
          y[j,i] = rnbinom( 1, size = 1/phi[j,cl[i]], mu = size_factors[i]*2^(b[j,cl[i]] + batch_eff[i]))
        } else{
          y[j,i] = rnbinom( 1, size = 1/phi[j,cl[i]], mu = size_factors[i]*2^(b[j,cl[i]]))
        }
      }
    }
  }
  result<-list(y=y,z=z,batch_g=batch_g)
  return(result)
}

# Simulate data with gene-wise dispersion parameters
simulate_data_g=function(n,k,g,init_pi,b,size_factors,distrib,phi=rep(0,times=g),
                         batch_coef=1,batch=rep(0,n),batch_g=NA){     # batch of 0: no batch effect, batch of 1: yes effect
  batch_eff = batch_coef*batch
  y<-matrix(rep(0,times=g*n),nrow=g)
  z = rmultinom(n,1,init_pi)
  if(all(is.na(batch_g))){
    batch_g = sample(1:g, 0.5*g)   # Apply batch on 50% of genes
  }
  
  if(ncol(b)!=k){
    warning("Wrong order selected. Simulating based on correct order")
    k=ncol(b)
  }
  cl = rep(0,n)
  for(c in 1:k){
    cl[z[c,]==1] = c
  }
  if(distrib=="poisson"){
    for(j in 1:g){
      for(i in 1:n){
        if(j %in% batch_g){
          y[j,i] = rpois( 1, lambda = size_factors[i]*2^(b[j,cl[i]] + batch_eff[i]))
        }else{
          y[j,i] = rpois( 1, lambda = size_factors[i]*2^(b[j,cl[i]]))
        }
      }
    }
  } else if(distrib=="nb"){
    for(j in 1:g){
      for(i in 1:n){
        if(j %in% batch_g){
          y[j,i] = rnbinom( 1, size = 1/phi[j], mu = size_factors[i]*2^(b[j,cl[i]] + batch_eff[i]))
        } else{
          y[j,i] = rnbinom( 1, size = 1/phi[j], mu = size_factors[i]*2^(b[j,cl[i]]))
        }
      }
    }
  }
  result<-list(y=y,z=z,batch_g=batch_g)
  return(result)
}

NB.GOF = function(y,size_factors=rep(1,ncol(y)),nsim=1000){
  # Algorithm based on Mi et al 2015
  n = ncol(y)
  g = nrow(y)
  R = nsim
  r0 = rep(0,times=n)
  pMC = rep(0,times=g)
  
  for(j in 1:g){
    start = Sys.time()
    cat("gene",j,"/",g,"\n")
    fit0 = glm.nb(as.numeric(y[j,]) ~ 1 + offset(log2(size_factors)),trace=0)
    #fit0 = glm.nb(y[j,] ~ 1)
    
    r0 = residuals(fit0,type="pearson")
    theta0 = fit0$theta
    coef0 = fit0$coefficients
    mu0 = 2^(coef0)
    rh = matrix(0,nrow=R,ncol=n)
    
    for(h in 1:R){
      yh = rnbinom(n,size=theta0,mu=mu0)
      #rh[h,] = (yh-mu0)/sqrt(mu0+mu0^2/theta0)            # referring to model fit0
      fith = glm.nb(yh ~ 1)      # Fitting a new model make sense??
      rh[h,] = residuals(fith,type="pearson")
      rh[h,]=rh[h,order(rh[h,])]
    }
    
    r2.5 = rep(0,n)
    r50 = rep(0,n)
    r97.5 = rep(0,n)
    
    for(i in 1:n){
      r2.5[i] = quantile(rh[,i],.025)
      r50[i] = quantile(rh[,i],.5)
      r97.5[i] = quantile(rh[,i],.975)
    }
    
    r0=r0[order(r0)]
    # plot(r0,r50,ylim=c(-3,4))
    # points(r0,r2.5,col="red")
    # points(r0,r97.5,col="blue")
    
    d0 = sum((r0-r50)^2)
    dh = rep(0,R)
    for(h in 1:R){
      dh[h] = sum((rh[h,]-r50)^2)
    }
    
    pMC[j] = (sum((dh >= d0)^2)+1) / (R+1)
    cat("pval =",pMC[j],"\n")
    end = Sys.time()
    cat("time_elap =",as.numeric(end-start),"\n")
  }
  return(pMC)
}

sim.iCluster = function(y,true_clusters,
                        ncores=10,n.lambda=25,K_min=2,K_max=7      # K_search = 2:7 for simulations, 2:15 for real data
                        ){
  
  # list of K to search over
  K_search=c(K_min:K_max)
  
  # iClusterPlus #
  iClust_OS <- list()
  for(c in (K_search-1)){
    cv.fit = tune.iClusterPlus(cpus=ncores,dt1=t(y),type="poisson",K=c,n.lambda=n.lambda,maxiter=20)
    iClust_OS[[c]] = cv.fit
  }
  BIC_mat = getBIC(iClust_OS)
  dev_mat = getDevR(iClust_OS)
  
  lambda_ids = apply(BIC_mat,2,which.min)   # indices of optimal lambda for each number of clusters
  devs=rep(0,length(K_search))
  for(c in (K_search-1)){
    devs[c] = dev_mat[lambda_ids[c],c]      # Deviance at optimal lambda for each cluster
  }                                         # Order selection is done by plateauing deviance value
     
  lambda_vals = iClust_OS[[1]]$lambda[lambda_ids]
  
  #max_k = which.max(devs)                # Problem: deviance always increases for data with any amount of noise
  
  dev_inc = rep(0,length(K_search)-1)
  for(i in 1:length(dev_inc)){
    dev_inc[i] = (devs[i+1]-devs[i])/devs[i]         # percent increase of POD
  }
  max_k=which(dev_inc<0.05)[1]                       # optimal cluster selected at index right before (less than 5% increase in POD) is observed
  if(is.na(max_k)){max_k=K_search[length(K_search)-1]}
  max_lambda = lambda_vals[max_k]
  
  iClust_fit <- iClusterPlus(dt1=t(y),type="poisson",lambda=max_lambda,K=max_k,maxiter=10)
  
  ARI = adjustedRandIndex(true_clusters,iClust_fit$clusters)
  
  results<-list(K=max_k+1,
                lambda=max_lambda,
                ARI=ARI,
                cls=iClust_fit$clusters,
                iClust_fit = iClust_fit)
  
  return(results)
}


sim.predict <- function(X,fit,new_dat,new_SF,true_clusters){
  res = predictions(X=X,fit=fit, newdata=new_dat, new_sizefactors=new_SF)
  pred_acc=adjustedRandIndex(res$final_clusters,true_clusters)
  return(list(pred_acc=pred_acc,sim_cls=res$final_clusters))
}

sim.NB.MClust = function(y,K_search){
  fit=NB.MClust(Count=t(y),K=K_search)
  
  max_k=fit$K
  cls=fit$cluster
  coefs = fit$parameters$mu
  phi = 1/fit$parameters$theta
  PPs = fit$parameters$posterior
  pi = fit$parameters$prior
  
  results = list(fit=fit,
                 max_k=max_k,
                 cls=cls,
                 coefs=coefs,
                 phi=phi,
                 PPs=PPs,
                 pi=pi)
  return(results)
}

# Function to perform EM on simulated data
sim.EM<-function(true.K, fold.change, num.disc, g, n, 
                 distrib,method="EM",filt_quant = 0.2,filt_method=c("pval","mad","none"),
                 sim_disp="gene",disp="gene",fixed_parms=F, fixed_coef=8,fixed_phi=0.35,
                 ncores=10,nsims=ncores,iCluster_compare=T,penalty=T,
                 sim_batch=F, p_batch=0.5, batch_coef=fold.change,
                 adj_batch=F){
  
  # disp: "gene" or "cluster"
  # low: coef 3.75-3.84, phi 0.13-0.15
  # med: coef 6.59-6.62, phi 0.32-0.38
  # high: coef 7.84-7.85, phi 1.00-1.32
  
  # Fixed phi: scalar for gene-wise, vector of length K for cluster-wise
  
  sim = nsims       # number of sims (set eq to number of cores for now)
  
  if(!fixed_parms){
    dir_name = sprintf("Sim_%d_%d_%d_%f_%f_%s_disp%s_batch%f_sim%s_adj%s",
                       true.K,n,g,fold.change,num.disc,distrib,sim_disp,p_batch,sim_batch,adj_batch)
  } else{
    if(length(fixed_phi)==1){
      dir_name = sprintf("Sim_%d_%d_%d_%f_%f_%s_fixed_disp%s_%f_%f_batch%f_sim%s_adj%s",
                         true.K,n,g,fold.change,num.disc,distrib,sim_disp,fixed_coef,fixed_phi,p_batch,sim_batch,adj_batch)
    }else{
      dir_name = sprintf("Sim_%d_%d_%d_%f_%f_%s_fixed_disp%s_%f_%s_batch%f_sim%s_adj%s",
                         true.K,n,g,fold.change,num.disc,distrib,sim_disp,fixed_coef,paste(fixed_phi,collapse="_"),p_batch,sim_batch,adj_batch)
    }
  }
  ifelse(!dir.exists(sprintf("Diagnostics/%s",dir_name)),
         dir.create(sprintf("Diagnostics/%s",dir_name)),
         FALSE)
  
  # max n = 100, max #
  
  # if(distrib=="poisson"){
  #   source("Pan EM.R")
  # } else if(distrib=="nb"){
  #   source("NB Pan EM par.R")
  # } else{
  #   print("no distrib input. Defaulting to Poisson")
  #   source("Pan EM.R")
  # }
  source("FSCseq.R")
  true_clusters<-NA        # TRUE clusters not known for real data
  
  fit_DESeq = fit_DESeq_intercept(init_y[1:g,1:n])
  init_size_factors <- fit_DESeq$size_factors
  init_norm_y <- fit_DESeq$norm_y
  
  # Unpenalized run to find initial cluster estimates based on K=k
  k=true.K
  if(!fixed_parms){
    res_init<-FSCseq(y=init_y,k=k,lambda=0,alpha=0,size_factors=init_size_factors,norm_y=init_norm_y,true_clusters=true_clusters,prefix="init",dir=dir_name,method=method,disp=sim_disp)
    init_coefs<-res_init$coefs              # save init estimates for coefs & pi
    init_phi<-res_init$phi
  } else{
    # fixed coefs and phi
    init_coefs <- matrix(fixed_coef,nrow=g,ncol=k)
    if(sim_disp=="gene"){
      init_phi <- rep(fixed_phi,g)
    } else{ init_phi <- matrix(fixed_phi,nrow=g,ncol=k,byrow=T) }
  }
  
  sim_coefs<-matrix(rep(rowSums(init_coefs)/k,times=k),ncol=k)
  fold.change<-fold.change
  nondisc.fold.change<-0         # fixed nondisc fold change
  tt<-floor(num.disc*g)
  
  down_disc = sample(1:tt,floor(tt/2),replace=F)          # Randomly apply LFC up on half of the genes, and down on half of the genes
  LFC_apply_mat =matrix(0,nrow=tt,ncol=k)
  LFC_apply = fold.change*(c(0:(k-1))+rep((1-k)/2,times=k))
  for(t in 1:tt){
    if(t %in% down_disc){
      LFC_apply_mat[t,]=-LFC_apply
    } else{
      LFC_apply_mat[t,]=LFC_apply
    }
  }
  sim_coefs[1:tt,]=sim_coefs[1:tt,] + LFC_apply_mat
  sim_pi<-rep(1/true.K,times=true.K)
  
  
  
  sink(file=sprintf("Diagnostics/%s/sim_parms_%s_%s.txt",dir_name,method,disp))
  cat("SIMULATED CLUSTER PROPORTIONS:\n")
  cat(sim_pi)
  cat("\n==========================================")
  cat("SIZE FACTORS:\n")
  cat(init_size_factors)
  cat("\n==========================================")
  cat("SIMULATED COEFFICIENTS:\n")
  write.table(sim_coefs,quote=F)
  cat("\n==========================================")
  cat("SIMULATED DISPERSION PARMS:\n")
  write.table(init_phi,quote=F)
  cat("\n==========================================")
  sink()
  
  
  #### SIMULATIONS ####
  
  # Simulations to find K (Order Selection)
  
  all_data <- list(list())
  
  for(ii in 1:sim){
    
    # Simulate data based on initial estimates/estimate size factors
    ## to simulate phi to be very small (fixed +10 extrapoisson variation)
    #init_phi=10/exp(rowMeans(sim_coefs))^2
    #init_phi = rchisq(g,2)
    #init_phi = rep(0,g)
    
    if(sim_batch){
      set.seed(ii)
      batch = apply(rmultinom(n,1,c(1-p_batch,p_batch)),2,which.max) - 1
    } else{
      batch = rep(0,n)
    }
    
    set.seed(2*ii)
    if(sim_disp=="cluster"){    # check for whether init_phi is of dimension 1
      sim.dat<-simulate_data(n=n,k=true.K,g=g,init_pi=sim_pi,b=sim_coefs,size_factors=init_size_factors,distrib=distrib,phi=init_phi,
                             batch_coef=batch_coef,batch=batch) # cluster-wise disp param
    } else{
      sim.dat<-simulate_data_g(n=n,k=true.K,g=g,init_pi=sim_pi,b=sim_coefs,size_factors=init_size_factors,distrib=distrib,phi=init_phi,
                               batch_coef=batch_coef,batch=batch) # gene-specific disp param
    }
    y<-sim.dat$y
    z<-sim.dat$z
    batch_g = rep(FALSE,g)
    batch_g[sim.dat$batch_g] <- TRUE
    true_clusters<-rep(0,times=n)
    for(i in 1:n){
      true_clusters[i]<-which(z[,i]==1)
    }
    norm_y = y
    for(i in 1:n){
      norm_y[,i] = y[,i]/init_size_factors[i]
    }
    true_disc=c(rep(TRUE,tt),rep(FALSE,(g-tt)))
    
    # # Filtering
    # idx <- rowMedians(norm_y) >= 20
    # y <- y[idx,]
    # norm_y <- norm_y[idx,]
    # true_disc <- true_disc[idx]
    
    # No filtering
    filt_ids = rep(T,g) # No filtering
    
    if(filt_method=="pval"){
      pvals = NB.GOF(y=y,size_factors=size_factors,nsim=1000)
      #FDR_pvals = p.adjust(pvals,"fdr")
      # pre-filtering by pval
      #filt_ids = (pvals <= pval_thresh)
      filt_ids = pvals <= quantile(pvals,filt_quant)
    } else if(filt_method=="mad"){
      mads = rep(0,g)
      for(j in 1:g){
        mads[j] = mad(log(norm_y[j,]+0.1))
      }
      filt_ids = mads >= quantile(mads,1-filt_quant)
    }
    idx = filt_ids
    
    y=y[idx,]
    norm_y=norm_y[idx,]
    true_disc=true_disc[idx]
    batch_g=batch_g[idx]
    
    all_data[[ii]]<-list(y=y,
                         true_clusters=true_clusters,
                         size_factors=init_size_factors,
                         norm_y=norm_y,
                         true_disc=true_disc, batch=batch, batch_g=batch_g
                         ,gene_id=idx
                        )
  }
    
  cat(paste(sim," datasets simulated \n"))
  
  # Function to run simulation in parallel
  sim.run = function(ii){
    
    y = all_data[[ii]]$y
    true_clusters = all_data[[ii]]$true_clusters
    true_disc = all_data[[ii]]$true_disc
    idx = all_data[[ii]]$gene_id
    
    filt_sens = sum(idx[1:tt])/tt
    filt_falsepos = sum(idx[(tt+1):g])/(g-tt)
    
    true_norm_y = all_data[[ii]]$norm_y               # ACTUAL norm_y based on simulated size factors
    true_size_factors = all_data[[ii]]$size_factors
    
    batch=all_data[[ii]]$batch
    batch_g=all_data[[ii]]$batch_g
    
    if(adj_batch){
      X=matrix(batch,ncol=1)
    }else{
      X=NA
    }
    
    fit=fit_DESeq_intercept(y,calc_vsd=F,calc_rld=F)       # ESTIMATE based on DESeq2
    size_factors=fit$size_factors
    norm_y=fit$norm_y
    rm(fit)
    
    sink(sprintf("Diagnostics/%s/progress%d_%s_%s.txt",dir_name,ii,method,disp))
    # Order selection
    K_search=c(1:7)
    list_BIC=matrix(0,nrow=length(K_search),ncol=2)
    list_BIC[,1]=K_search
    print(paste("Dataset",ii,"Order Selection:"))
    list_res = list()
    
    for(aa in K_search){
      pref = sprintf("order%d",ii)
      start=as.numeric(Sys.time())
      res<-FSCseq(X=X,y=y,k=list_BIC[aa,1],lambda=0,alpha=0,size_factors=size_factors,norm_y=norm_y,
            true_clusters=true_clusters,true_disc=true_disc,prefix=pref,dir=dir_name,method=method,disp=disp)  # alpha = 0: all L1. No penalty here
      end=as.numeric(Sys.time())
      list_BIC[aa,2]<-res$BIC
      # if(list_BIC[aa,1]==true.K){
      #   compare_res = res
      # }
      print(list_BIC[aa,])
      print(paste("Time:",end-start,"seconds"))
      list_res[[aa]]=res
    }
    max_k=list_BIC[which.min(list_BIC[,2]),1]
    unpen_BIC = min(list_BIC[,2])
    
    # sink(file=sprintf("Diagnostics/%s/%s_%s_final%d_order.txt",dir_name,method,disp,ii))
    # cat(paste("True order:",true.K,"\n"))
    # cat(paste("Optimal order selected:",max_k,"\n"))
    # cat("RUN WITH CORRECT ORDER:\n")
    # MSE_coefs = sum((compare_res$coefs[,order(compare_res$coefs[1,])] - sim_coefs[idx,])^2)/(sum(idx)*true.K)
    # if(is.null(ncol(init_phi))){
    #   MSE_phi = sum((init_phi[idx]-compare_res$phi)^2)/(sum(idx)*true.K) # test
    # } else{
    #   MSE_phi = sum((init_phi[idx,]-compare_res$phi)^2)/(sum(idx)*true.K) # test
    # }
    # cat(paste("ARI:",adjustedRandIndex(compare_res$final_clusters,true_clusters),"\n"))
    # cat(paste("MSE of true vs discovered coefs:",MSE_coefs,"\n"))
    # cat(paste("MSE of true vs discovered phi:",MSE_phi,"\n"))
    # cat(paste("% of correctly ID'ed disc genes:",sum(!compare_res$nondiscriminatory==true_disc)/sum(true_disc),"\n"))
    # cat(paste("PPs (n x k):\n"))
    # write.table(t(compare_res$wts),quote=F,col.names=F)
    # 
    # sink()
    # 
    # pdf(file=sprintf("Diagnostics/%s/%s_%s_final%d_order.pdf",dir_name,method,disp,ii))
    # for(c in 1:true.K){
    #   cl_ids = true_clusters==c
    #   for(cc in 1:true.K){
    #     boxplot(compare_res$wts[cc,cl_ids],main=sprintf("Boxplot of PP for subjects of true cl%d being in cl%d",c,cc))
    #   }
    # }
    # annotation_col = data.frame(cbind(true_clusters,compare_res$final_clusters))
    # colnames(annotation_col)=c("True","Derived")
    # rownames(annotation_col)=c(1:ncol(norm_y))
    # colnames(norm_y)=c(1:ncol(norm_y))
    # annotation_col2 = annotation_col[order(true_clusters),]
    # pheatmap(log(norm_y[,order(compare_res$final_clusters)]+0.1),cluster_cols = F,scale="row",annotation_col = annotation_col2)
    # dev.off()
    
    # Grid search
    
    # # Use prefiltering just on the order selection step? This will reset all genes
    # y = all_data[[ii]]$y
    # true_clusters = all_data[[ii]]$true_clusters
    # norm_y = all_data[[ii]]$norm_y
    # true_disc = all_data[[ii]]$true_disc
    # idx = rep(T,nrow(y))
    
    if(penalty){
      print(paste("Dataset",ii,"Grid Search:"))    # track iteration
      
      #create matrix for grid search values
      lambda_search=seq(0.05,1,0.05)
      alpha_search=seq(0,0.5,0.05)        # alpha can't = 1
      
      list_BIC=matrix(0,nrow=length(lambda_search)*length(alpha_search),ncol=3) # matrix of BIC's: one for each combination of penalty params 
      
      list_BIC[,1]=rep(lambda_search,each=length(alpha_search))
      list_BIC[,2]=rep(alpha_search,times=length(lambda_search))
      
      #list_pen_res = list()
      
      #search for optimal penalty parameters
      ##################################################### figure out why higher penalty errors out (lambda=.75, alpha=.5) but lower is fine (lambda=.75, alpha=.45)####
      for(aa in 1:nrow(list_BIC)){
        pref = sprintf("grid%d",ii)
        start = as.numeric(Sys.time())
        res<-FSCseq(X=X,y=y,k=max_k,lambda=list_BIC[aa,1],alpha=list_BIC[aa,2],size_factors=size_factors,norm_y=norm_y,
              true_clusters=true_clusters,true_disc=true_disc,disp=disp,prefix=pref,dir=dir_name,method=method,
              init_cls=list_res[[max_k]]$final_clusters,init_parms=T,init_coefs=list_res[[max_k]]$coefs,init_phi=list_res[[max_k]]$phi)
        end = as.numeric(Sys.time())
        list_BIC[aa,3]<-res$BIC
        print(list_BIC[aa,])
        print(paste("Time:",end-start,"seconds"))
        #list_pen_res[[aa]]=res
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
      
      if(pen_BIC>=unpen_BIC){
        max_lambda=0
        max_alpha=0
        print("Algorithm selected unpenalized model")
      } else{
        print(paste("Penalized model with lambda=",max_lambda,"and alpha=",max_alpha))
      }
    } else {
      max_lambda=0
      max_alpha=0
      print(paste("No Penalization"))
    }
    
    # Final run with optimal parameters
    pref = sprintf("final%d",ii)
    start = as.numeric(Sys.time())
    res<-FSCseq(X=X,y=y,k=max_k,lambda=max_lambda,alpha=max_alpha,size_factors=size_factors,norm_y=norm_y,
          true_clusters=true_clusters,true_disc=true_disc,prefix=pref,dir=dir_name,method=method,disp=disp,
          init_cls=list_res[[max_k]]$final_clusters,init_parms=T,init_coefs=list_res[[max_k]]$coefs,init_phi=list_res[[max_k]]$phi)
    end = as.numeric(Sys.time())
    res$time_elap = end-start
    
    print("Optimal EM run complete")
    
    cls_EM=res$final_clusters
    
    
    # iCluster+
    if(iCluster_compare){
      n.lambda=25
      iCluster_res = sim.iCluster(round(norm_y,0),true_clusters,ncores=1,n.lambda=n.lambda,K_min=2,K_max=7)
      # if(iCluster_res$K != true.K){
      #   cv.fit = tune.iClusterPlus(cpus=1,dt1=t(y),type="poisson",K=true.K-1,alpha=1,n.lambda=n.lambda,scale.lambda=1,maxiter=20)
      #   BIC_mat = matrix(0,nrow=25,ncol=2)
      #   for(i in 1:n.lambda){
      #     BIC_mat[i,1] = cv.fit$lambda[i]
      #     BIC_mat[i,2] = cv.fit$fit[[i]]$BIC
      #   }
      #   max_lambda = BIC_mat[which.min(BIC_mat[,2]),1]
      #   
      #   iClust_fit <- iClusterPlus(dt1=t(y),type="poisson",lambda=max_lambda,alpha=1,K=true.K-1,maxiter=10)
      # }
      K_iClust = iCluster_res$K
      cls_iClust = iCluster_res$cls
      ARI_iClust = iCluster_res$ARI
    } else{iCluster_res=NA}
    rm("iCluster_res")
    
    print("iCluster run complete")
    
    
    # Predictions:
    y_pred=NA
    pred_acc=NA
    true_clusters_pred=NA
      res_pred = res
    
      #n_pred = floor(0.1*n)           # simulate data for 10% of original n
      n_pred = 20                     # FIXED number of samples to predict
      SF_pred = size_factors[1:n_pred]
      
      if(sim_batch){
        set.seed(ii)
        batch_pred = apply(rmultinom(n_pred,1,c(1-p_batch,p_batch)),2,which.max) - 1
      } else{
        batch_pred = rep(0,n_pred)
      }
      
      if(adj_batch){
        X_pred=matrix(batch_pred,ncol=1)
      }else{
        X_pred=NA
      }
      
      if(sim_disp=="cluster"){    # check for whether init_phi is of dimension 1
        sim.dat<-simulate_data(n=n_pred,k=true.K,g=g,init_pi=sim_pi,b=sim_coefs,size_factors=SF_pred,distrib=distrib,phi=init_phi,
                               batch_coef=batch_coef,batch=batch_pred,batch_g=batch_g) # cluster-wise disp param
      } else{
        sim.dat<-simulate_data_g(n=n_pred,k=true.K,g=g,init_pi=sim_pi,b=sim_coefs,size_factors=SF_pred,distrib=distrib,phi=init_phi,
                                 batch_coef=batch_coef,batch=batch_pred,batch_g=batch_g) # gene-specific disp param
      }
      y_pred<-sim.dat$y
      z_pred<-sim.dat$z
      true_clusters_pred<-rep(0,times=n_pred)
      
      # same pre-filtering for prediction dataset:
      y_pred = y_pred[idx,]
      
      for(i in 1:n_pred){
        true_clusters_pred[i]<-which(z_pred[,i]==1)
      }
      
      # # match cluster ID'sbased on SSE's of true vs estimated coefficients
      # SSEs=matrix(0,nrow=true.K,ncol=true.K)
      # subs_sim_coefs=sim_coefs[idx,]
      # for(c in 1:true.K){
      #   for(cc in 1:true.K){
      #     SSEs[c,cc] = sum(abs(subs_sim_coefs[,c]-res_pred$coefs[,cc]))
      #   }
      # }
      # new_cl_ids = rep(0,true.K)
      # for(c in 1:true.K){
      #   new_cl_ids[c]=which.min(SSEs[c,])
      # }
      # true_clusters_pred=new_cl_ids[true_clusters_pred]
      
      
      fit= sim.predict(X_pred,res_pred,y_pred,SF_pred,true_clusters=true_clusters_pred)
      pred_acc=fit$pred_acc
    
    print("Prediction complete")
    
    set.seed(123)
    results_hc=list()
    selected <- c( "kl", "ch", "hartigan",  "cindex", "db", "silhouette", "duda", "pseudot2", "beale", "ratkowsky", "ball", "ptbiserial", "gap", "mcclain", "gamma", "gplus", "tau", "dunn", "hubert", "sdindex", "dindex", "sdbw")
    max_K_hc_avail=0
    
    for (i in 1:length(selected)) {
      cat(paste("HC order select method: ",selected[i],"\n"))
      pdf(sprintf("Diagnostics/%s/hc_OS%d.pdf",dir_name,ii))
      results_hc[[i]] <- try(NbClust(t(log(norm_y+0.1)), min.nc=min(K_search), max.nc=max(K_search), method="average", index=selected[i]))
      dev.off()
      max_K_hc_avail=max_K_hc_avail+(length(results_hc[[i]])>1)
    }
    max_K_hcs = rep(NA,max_K_hc_avail)
    for(i in 1:max_K_hc_avail){
      if(length(results_hc[[i]])>1){
        max_K_hcs[i] = try(results_hc[[i]]$Best.nc[1])
      }
    }
    K_hc = as.numeric(unique(max_K_hcs)[order(unique(max_K_hcs))][which.max(table(max_K_hcs))])
    cat(paste("HC K: ",K_hc, "\n"))
    
    pam_fit = list()
    sil.vals = rep(NA,max(K_search))
    for(c in K_search){
      if(c==1) next         # can't search k=1
      pam_fit[[c]] = pam(t(log(norm_y+0.1)),c)
      sil.vals[c] = pam_fit[[c]]$silinfo$avg.width
      cat(paste("K-med order select avg silhouette: ", sil.vals[c],"\n"))
    }
    K_med = which.max(sil.vals)
    cat(paste("K-med K: ",K_med, "\n"))
    
    cls_hc = cutree(hclust(as.dist(1-cor(norm_y, method="spearman")),method="average"),K_hc)
    cat("cls_hc finished\n")
    cls_med = pam(t(log(norm_y+0.1)),K_med)$cluster
    cat("cls_km finished\n")
    
    
    # NB.MClust
    NBMB_fit=NA
    K_NBMB=0
    cls_NBMB=rep(0,n)
    tryCatch({
      NBMB_fit = sim.NB.MClust(round(norm_y,0),2:7)
      K_NBMB = NBMB_fit$max_k
      cls_NBMB =NBMB_fit$cls
    
      cat(paste("NB.MClust K:",K_NBMB,"\n"))
      cat("Finished NB.MClust\n")
    },error=function(e){
      cat("NB.MClust errored out (glm.nb fit did not converge). Excluding from analysis\n")
    })
    
    # Mclust (log, variance-stabilizing, and rlog transforms)
    fit = fit_DESeq_intercept(y,calc_vsd=T,calc_rld=T)
    log_norm_y=log(norm_y+0.1)
    if(!is.na(fit$dds)){
      vsd_y = assay(fit$vsd)
      rld_y = assay(fit$rld)
    }
    cat("Finished performing vsd/rld transformations\n")
    
    K_vsd_mclust=NA
    K_rld_mclust=NA
    cls_vsd_mclust=NA
    cls_rld_mclust=NA
    
    mclust_log_fit=Mclust(t(log_norm_y),G=1:7)
    K_log_mclust = mclust_log_fit$G
    cls_log_mclust = mclust_log_fit$classification
    rm("mclust_log_fit")
    cat("mclust log transform order: ",K_log_mclust,"\n")
    
    if(!is.na(fit$dds)){
      mclust_vsd_fit=Mclust(t(vsd_y),G=1:7)
      K_vsd_mclust = mclust_vsd_fit$G
      cls_vsd_mclust = mclust_vsd_fit$classification
      rm("mclust_vsd_fit")
      cat("mclust vsd transform order: ",K_vsd_mclust,"\n")
      
      mclust_rld_fit=Mclust(t(rld_y),G=1:7)
      K_rld_mclust = mclust_rld_fit$G
      cls_rld_mclust = mclust_rld_fit$classification
      rm("mclust_rld_fit")
      cat("mclust rld transform order: ",K_rld_mclust,"\n")
    }
    
    rm("fit")
    
    
    print(paste("Time:",res$time_elap,"seconds"))
    print(paste("Dataset ",ii,"complete"))
    
    sink()
    
    all_Ks = c(max_k,K_iClust,K_NBMB,K_log_mclust,K_vsd_mclust,K_rld_mclust,K_hc,K_med)
    names(all_Ks) = c("EM","iClust","NBMB","logMC","vsdMC","rldMC","HC","KM")
    all_cls = cbind(cls_EM,cls_iClust,cls_NBMB,
          cls_log_mclust,cls_vsd_mclust,cls_rld_mclust,
          cls_hc,cls_med)
    colnames(all_cls) = c("EM","iClust","NBMB","logMC","vsdMC","rldMC","HC","KM")
    
    print("all_Ks and all_cls created")
    
    results=list(res=res,
                 max_k=max_k,
                 max_lambda=max_lambda,
                 max_alpha=max_alpha,
                 true_clusters=true_clusters,
                 true_disc=true_disc,
                 pred_acc=pred_acc,
                 y_pred=y_pred,
                 true_clusters_pred=true_clusters_pred,
                 all_Ks=all_Ks,all_cls=all_cls,
                 filt_sens=filt_sens,
                 filt_falsepos=filt_falsepos, batch=batch, batch_pred=batch_pred)
    return(results)
  }

  ## ADD PARALLELIZATION HERE ##
  cl<-makeCluster(ncores)
  clusterExport(cl=cl,varlist=c(ls(),"fit_DESeq_intercept","FSCseq","EM_run","logsumexpc","SCAD_soft_thresholding","lasso_soft_thresholding","NB.GOF",
                                "simulate_data","simulate_data_g","sim.iCluster","sim.predict","predictions","sim.NB.MClust"),envir=environment())
  clusterEvalQ(cl,{
    library(stats)
    library(MASS)
    library(permute)
    library(Rcpp)
    library(RcppArmadillo)
    library(mclust)
    library(pryr)
    library(pheatmap)
    library(DESeq2)
    library(iClusterPlus)
    library(cluster)
    library(NbClust)
    library(NB.MClust)
    library(Biobase)
    sourceCpp("M_step.cpp")
    })
  
  # Store all results in list par_sim_res
  par_sim_res=list(list())
  par_sim_res<-parLapply(cl, 1:sim, sim.run)
  stopCluster(cl)
  # # non parallel computation (for debugging)
  # for(ii in 1:sim){
  #   par_sim_res[[ii]] = sim.run(ii)
  # }
  
  print("Finished parallel computations")
  
  
  all_sim_data <- list(list())
  
  # Initialize EM results
  temp_ks<-rep(0,times=sim)
  temp_lambdas<-rep(0,times=sim)
  temp_alphas<-rep(0,times=sim)
  #temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
  #temp_coefs<-list()
  temp_disc<-rep(0,times=sim)
  temp_ARI<-rep(0,times=sim)
  temp_sensitivity<-rep(0,times=sim)
  temp_falsepos<-rep(0,times=sim)
  temp_pred_acc<-rep(0,times=sim)
  
  # Initialize all other results Ks
  temp_K_iClust = rep(0,sim)
  temp_K_hc = rep(0,sim)
  temp_K_med = rep(0,sim)
  temp_K_NBMB = rep(0,sim)
  temp_K_log_mclust = rep(0,sim)
  temp_K_vsd_mclust = rep(0,sim)
  temp_K_rld_mclust = rep(0,sim)
  
  # Initialize all other results ARIs
  temp_ARI_iClust = rep(0,sim)
  temp_ARI_hc = rep(0,sim)
  temp_ARI_med = rep(0,sim)
  temp_ARI_NBMB = rep(0,sim)
  temp_ARI_log_mclust = rep(0,sim)
  temp_ARI_vsd_mclust = rep(0,sim)
  temp_ARI_rld_mclust = rep(0,sim)
  
  
  # temp_sil_hc = rep(0,sim)
  # temp_sil_med = rep(0,sim)
  # temp_sil_EM = rep(0,sim)
  # temp_sil_iClust = rep(0,sim)

  temp_filt_sens = rep(0,sim)
  temp_filt_falsepos = rep(0,sim)
  
  all_res = list()
  all_Ks = list()
  all_cls = list()
  all_batch = list()
  all_batch_pred = list()

  # Summarize results
  for(ii in 1:sim){
    true_clusters=par_sim_res[[ii]]$true_clusters
    true_disc=par_sim_res[[ii]]$true_disc
    all_Ks[[ii]]=par_sim_res[[ii]]$all_Ks
    all_cls[[ii]]=par_sim_res[[ii]]$all_cls
    
    all_batch[[ii]]=par_sim_res[[ii]]$batch
    all_batch_pred[[ii]]=par_sim_res[[ii]]$batch_pred
    
    # Store EM results
    res=par_sim_res[[ii]]$res
    all_res[[ii]]=res
    temp_ks[ii]=all_Ks[[ii]]["EM"]
    temp_lambdas[ii] = par_sim_res[[ii]]$max_lambda
    temp_alphas[ii] = par_sim_res[[ii]]$max_alpha
    temp_disc[ii]<-sum(!res$nondiscriminatory)
    temp_ARI[ii]<-adjustedRandIndex(true_clusters,res$final_clusters)
    temp_pred_acc[ii] <- par_sim_res[[ii]]$pred_acc
    if(tt>0){
      temp_sensitivity[ii]<-sum(!res$nondiscriminatory[true_disc])/tt                   # tt = # of disc genes simulated (floor(num.disc*g))
    } else {temp_sensitivity[ii]<-NA}
    if(tt<g){
      temp_falsepos[ii]<-sum(!res$nondiscriminatory[!true_disc])/(g-tt)
    } else {temp_falsepos[ii]<-NA}         # take into account genes omitted b/c rowSum of count was <100 (such genes are assumed nondiscriminatory)
    temp_filt_sens[ii] = par_sim_res[[ii]]$filt_sens
    temp_filt_falsepos[ii] = par_sim_res[[ii]]$filt_falsepos
    all_sim_data[[ii]] = list(par_sim_res[[ii]]$y_pred,par_sim_res[[ii]]$true_clusters_pred)
    
    # Store all other Ks
    temp_K_iClust[ii] = all_Ks[[ii]]["iClust"]
    temp_K_hc[ii] = all_Ks[[ii]]["HC"]
    temp_K_med[ii] = all_Ks[[ii]]["KM"]
    temp_K_NBMB[ii] = all_Ks[[ii]]["NBMB"]
    temp_K_log_mclust[ii] = all_Ks[[ii]]["logMC"]
    temp_K_vsd_mclust[ii] = all_Ks[[ii]]["vsdMC"]
    temp_K_rld_mclust[ii] = all_Ks[[ii]]["rldMC"]
    
    temp_ARI_iClust[ii] = adjustedRandIndex(all_cls[[ii]][,"iClust"],true_clusters)
    temp_ARI_hc[ii] = adjustedRandIndex(all_cls[[ii]][,"HC"],true_clusters)
    temp_ARI_med[ii] = adjustedRandIndex(all_cls[[ii]][,"KM"],true_clusters)
    temp_ARI_NBMB[ii] = adjustedRandIndex(all_cls[[ii]][,"NBMB"],true_clusters)
    temp_ARI_log_mclust[ii] = adjustedRandIndex(all_cls[[ii]][,"logMC"],true_clusters)
    temp_ARI_vsd_mclust[ii] = adjustedRandIndex(all_cls[[ii]][,"vsdMC"],true_clusters)
    temp_ARI_rld_mclust[ii] = adjustedRandIndex(all_cls[[ii]][,"rldMC"],true_clusters)
    
  }
  
  #### AVERAGE EM RESULTS ####
  #mean_pi<-colSums(temp_pi)/sim
  #mean_coefs<-Reduce('+',temp_coefs)/sim
  mean_disc<-mean(temp_disc)/g
  mean_ARI<-mean(temp_ARI)
  mean_sensitivity<-mean(temp_sensitivity)
  mean_falsepos<-mean(temp_falsepos)
  K=mean(temp_ks)
  mean_pred_acc = mean(temp_pred_acc,na.rm=T)
  mean_order_pred_acc = mean(temp_pred_acc[temp_ks==true.K],na.rm=T)
  
  #### AVERAGE OTHER Ks ####
  K_iClust = mean(temp_K_iClust)         # range: 1:7
  K_HC = mean(temp_K_hc)                 # range: 1:7
  K_KM = mean(temp_K_med)                # range: 2:7
  K_NBMB = mean(temp_K_NBMB[temp_K_NBMB!=0])             # range: 2:7
  K_log_MC = mean(temp_K_log_mclust)     # range: 1:7
  K_vsd_MC = mean(temp_K_vsd_mclust)     # range: 1:7
  K_rld_MC = mean(temp_K_rld_mclust)     # range: 1:7
  
  #### AVERAGE OTHER ARIs ####
  ARI_iClust = mean(temp_ARI_iClust)
  ARI_HC = mean(temp_ARI_hc)
  ARI_KM = mean(temp_ARI_med)
  ARI_NBMB = mean(temp_ARI_NBMB[temp_ARI_NBMB!=0])
  ARI_log_MC = mean(temp_ARI_log_mclust)
  ARI_vsd_MC = mean(temp_ARI_vsd_mclust)
  ARI_rld_MC = mean(temp_ARI_rld_mclust)
  
  #### Pre-filtering Sens/FP average ####
  filt_sens = mean(temp_filt_sens)
  filt_falsepos = mean(temp_filt_falsepos)
  
  # #### Silhouette values comparisons
  # mean_sil_hc = mean(temp_sil_hc)
  # mean_sil_med = mean(temp_sil_med)
  # mean_sil_EM = mean(temp_sil_EM)
  # mean_sil_iClust = mean(temp_sil_iClust)
  
  # Results for tabulation:
  
  results<-list(all_data=all_data,
                all_sim_data=all_sim_data,
                
                all_res=all_res,
                K=K,
                all_k=temp_ks,
                all_lambda=temp_lambdas,
                all_alpha=temp_alphas,
                all_ARI=temp_ARI,
                ARI=mean_ARI,
                all_disc=temp_disc,
                disc=mean_disc,
                all_sens=temp_sensitivity,
                sens=mean_sensitivity,
                all_falsepos=temp_falsepos,
                falsepos=mean_falsepos,
            
                all_pred_acc=temp_pred_acc,
                mean_pred_acc=mean_pred_acc,
                
                all_Ks=all_Ks,all_cls=all_cls,
                
                K_iClust=K_iClust, K_HC=K_HC, K_KM=K_KM, K_NBMB=K_NBMB, K_log_MC=K_log_MC, K_vsd_MC=K_vsd_MC, K_rld_MC=K_rld_MC,
                ARI_iClust=ARI_iClust, ARI_HC=ARI_HC, ARI_KM=ARI_KM, ARI_NBMB=ARI_NBMB, ARI_log_MC=ARI_log_MC, ARI_vsd_MC=ARI_vsd_MC, ARI_rld_MC=ARI_rld_MC,
                
                filt_sens=filt_sens,filt_falsepos=filt_falsepos,
                all_filt_sens=temp_filt_sens,
                all_filt_falsepos=temp_filt_falsepos,
                
                all_K_iClust=temp_K_iClust,all_K_HC=temp_K_hc,all_K_KM=temp_K_med,all_K_NBMB=temp_K_NBMB,all_K_log_MC=temp_K_log_mclust,all_K_vsd_MC=temp_K_vsd_mclust,all_K_rld_MC=temp_K_rld_mclust,
                all_ARI_iClust=temp_ARI_iClust,all_ARI_HC=temp_ARI_hc,all_ARI_KM=temp_ARI_med,all_ARI_NBMB=temp_ARI_NBMB,all_ARI_log_MC=temp_ARI_log_mclust,all_ARI_vsd_MC=temp_ARI_vsd_mclust,all_ARI_rld_MC=temp_ARI_rld_mclust,
                
                mean_order_pred_acc=mean_order_pred_acc,
                all_batch=all_batch,all_batch_pred=all_batch_pred
                
                )
  
  return(results)
}
  
  
  
  ###################################
  ###################################
  ###################################
  
  
