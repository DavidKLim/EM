#setwd("C:/Users/David/Desktop/Research/EM")

# NSCLC
#init_y<-read.table("init_y.txt")
#init_size_factors<-as.numeric(read.table("init_size_factors.txt")[,1])
#init_norm_y<-read.table("init_norm_y.txt")
#init_y<-cbind(init_y,init_y,init_y,init_y,init_y,init_y,init_y,init_y,init_y,init_y,init_y,init_y)

# BRCA
load('BRCA_env.RData')
init_y = y

library("stats")
library("data.table")
library("DESeq2")
library("mclust")
library("parallel")
library("pheatmap")
library("iClusterPlus")
library("cluster")



# Function to simulate data
simulate_data=function(n,k,g,init_pi,b,size_factors,distrib,phi=matrix(0,nrow=g,ncol=k)){
  y<-matrix(rep(0,times=g*n),nrow=g)
  z = rmultinom(n,1,init_pi)
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
        y[j,i] = rpois( 1, lambda = size_factors[i]*exp(b[j,cl[i]]))
      }
    }
  } else if(distrib=="nb"){
    for(j in 1:g){
      for(i in 1:n){
        y[j,i] = rnbinom( 1, size = 1/phi[j,cl[i]], mu = size_factors[i]*exp(b[j,cl[i]]))
      }
    }
  }
  result<-list(y=y,z=z)
  return(result)
}

# Simulate data with gene-wise dispersion parameters
simulate_data_g=function(n,k,g,init_pi,b,size_factors,distrib,phi=rep(0,times=g)){
  y<-matrix(rep(0,times=g*n),nrow=g)
  z = rmultinom(n,1,init_pi)
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
        y[j,i] = rpois( 1, lambda = size_factors[i]*exp(b[j,cl[i]]))
      }
    }
  } else if(distrib=="nb"){
    for(j in 1:g){
      for(i in 1:n){
        y[j,i] = rnbinom( 1, size = 1/phi[j], mu = size_factors[i]*exp(b[j,cl[i]]))
      }
    }
  }
  result<-list(y=y,z=z)
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
    fit0 = glm.nb(as.numeric(y[j,]) ~ 1 + offset(log(size_factors)),trace=0)
    #fit0 = glm.nb(y[j,] ~ 1)
    
    r0 = residuals(fit0,type="pearson")
    theta0 = fit0$theta
    coef0 = fit0$coefficients
    mu0 = exp(coef0)
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
                        ncores=10,n.lambda=25){
  
  # list of K to search over
  K_search = c(2:7)
  
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


sim.predict <- function(X,new_dat,new_SF,true_clusters){
  X = predictions(X=X, newdata=new_dat, new_sizefactors=new_SF)
  pred_acc=mean(X$final_clusters==true_clusters)
  return(list(pred_acc=pred_acc,sim_cls=X$final_clusters))
}



# Function to perform EM on simulated data
sim.EM<-function(true.K, fold.change, num.disc, g, n, 
                 distrib,method="EM",filt_quant = 0.2,filt_method=c("pval","mad","none"),
                 disp="gene",fixed_parms=F, fixed_coef=6.5,fixed_phi=0.35,
                 ncores=10,nsims=ncores,iCluster_compare=F,penalty=T){
  
  # disp: "gene" or "cluster"
  # low: coef 3.75-3.84, phi 0.13-0.15
  # med: coef 6.59-6.62, phi 0.32-0.38
  # high: coef 7.84-7.85, phi 1.00-1.32
  
  # Fixed phi: scalar for gene-wise, vector of length K for cluster-wise
  
  sim = nsims       # number of sims (set eq to number of cores for now)
  
  if(!fixed_parms){
    dir_name = sprintf("Sim_%d_%d_%d_%f_%f_%s",true.K,n,g,fold.change,num.disc,distrib)
  } else{
    if(length(fixed_phi)==1){
      dir_name = sprintf("Sim_%d_%d_%d_%f_%f_%s_fixed_%f_%f",true.K,n,g,fold.change,num.disc,distrib,fixed_coef,fixed_phi)
    }else{
      dir_name = sprintf("Sim_%d_%d_%d_%f_%f_%s_fixed_%f_%s",true.K,n,g,fold.change,num.disc,distrib,fixed_coef,paste(fixed_phi,collapse="_"))
    }
  }
  dir.create(sprintf("Diagnostics/%s",dir_name))
  
  # max n = 100, max #
  
  if(distrib=="poisson"){
    source("Pan EM.R")
  } else if(distrib=="nb"){
    source("NB Pan EM par.R")
  } else{
    print("no distrib input. Defaulting to Poisson")
    source("Pan EM.R")
  }
  true_clusters<-NA        # TRUE clusters not known for real data
  
  
  init_y<-init_y[1:g,1:n]
  row_names<-paste("gene",seq(g))
  col_names<-paste("subj",seq(n))
  cts<-as.matrix(init_y)
  rownames(cts)<-row_names
  colnames(cts)<-col_names
  coldata<-data.frame(matrix(paste("cl",true_clusters,sep=""),nrow=n))
  rownames(coldata)<-colnames(cts)
  colnames(coldata)<-"cluster"
  dds<-DESeqDataSetFromMatrix(countData = cts,
                              colData = coldata,
                              design = ~ 1)
  dds<-DESeq(dds)
  init_size_factors<-sizeFactors(dds)
  init_norm_y<-counts(dds,normalized=TRUE)
  
  # Unpenalized run to find initial cluster estimates based on K=k
  k=true.K
  if(!fixed_parms){
    X_init<-EM(y=init_y,k=k,lambda1=0,lambda2=0,tau=0,size_factors=init_size_factors,norm_y=init_norm_y,true_clusters=true_clusters,prefix="init",dir=dir_name,method=method,disp=disp)
    init_coefs<-X_init$coefs              # save init estimates for coefs & pi
    init_phi<-X_init$phi
  } else{
    # fixed coefs and phi
    init_coefs <- matrix(fixed_coef,nrow=g,ncol=k)
    if(disp=="gene"){
      init_phi <- rep(fixed_phi,g)
    } else{ init_phi <- matrix(fixed_phi,nrow=g,ncol=k,byrow=T) }
  }
  
  
  size_factors<-init_size_factors    # use this for all simulations
  
  sim_coefs<-matrix(rep(rowSums(init_coefs)/k,times=k),ncol=k)
  fold_change<-fold.change
  nondisc_fold_change<-0         # fixed nondisc fold change
  tt<-floor(num.disc*g)
  sim_coefs[1:tt,]<-matrix(rep( fold_change*(c(0:(k-1))+rep((1-k)/2,times=k)) ,times=tt),nrow=tt,byrow=TRUE)+sim_coefs[1:tt,]
  #sim_coefs[(tt+1):g,]<-matrix(rep( nondisc_fold_change*(c(0:(k-1))+rep((1-k)/2,times=k)) ,times=(g-tt)),nrow=(g-tt),byrow=TRUE)+sim_coefs[(tt+1):g,]         # nondisc fold change = 0 so this doesn't get changed
  sim_pi<-rep(1/true.K,times=true.K)
  
  
  
  sink(file=sprintf("Diagnostics/%s/sim_parms_%s_%s.txt",dir_name,method,disp))
  cat("SIMULATED CLUSTER PROPORTIONS:\n")
  cat(sim_pi)
  cat("\n==========================================")
  cat("SIZE FACTORS:\n")
  cat(size_factors)
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
    
    if(disp=="cluster"){    # check for whether init_phi is of dimension 1
      sim.dat<-simulate_data(n=n,k=true.K,g=g,init_pi=sim_pi,b=sim_coefs,size_factors=size_factors,distrib=distrib,phi=init_phi) # cluster-wise disp param
    } else{
      sim.dat<-simulate_data_g(n=n,k=true.K,g=g,init_pi=sim_pi,b=sim_coefs,size_factors=size_factors,distrib=distrib,phi=init_phi) # gene-specific disp param
    }
    y<-sim.dat$y
    z<-sim.dat$z
    true_clusters<-rep(0,times=n)
    for(i in 1:n){
      true_clusters[i]<-which(z[,i]==1)
    }
    norm_y = y
    for(i in 1:n){
      norm_y[,i] = y[,i]/size_factors[i]
    }
    true_disc=c(rep(TRUE,tt),rep(FALSE,(g-tt)))
    
    # # Filtering
    # idx <- rowMedians(norm_y) >= 20
    # y <- y[idx,]
    # norm_y <- norm_y[idx,]
    # true_disc <- true_disc[idx]
    
    # No filtering
    #idx = rep(T,g)
    
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
    
    all_data[[ii]]<-list(y=y,
                         true_clusters=true_clusters,
                         size_factors=size_factors,
                         norm_y=norm_y,
                         true_disc=true_disc
                         ,gene_id=idx
                        )
  }
    
  cat(paste(sim," datasets simulated \n"))
  
  # Function to run simulation in parallel
  sim.run = function(ii){
    
    y = all_data[[ii]]$y
    true_clusters = all_data[[ii]]$true_clusters
    norm_y = all_data[[ii]]$norm_y
    true_disc = all_data[[ii]]$true_disc
    idx = all_data[[ii]]$gene_id
    
    sink(sprintf("Diagnostics/%s/progress%d.txt",dir_name,ii))
    # Order selection
    K_search=c(2:7)
    list_BIC=matrix(0,nrow=length(K_search),ncol=2)
    list_BIC[,1]=K_search
    print(paste("Dataset",ii,"Order Selection:"))
    
    for(aa in 1:nrow(list_BIC)){
      pref = sprintf("order%d",ii)
      start=as.numeric(Sys.time())
      X<-EM(y=y,k=list_BIC[aa,1],lambda1=0,lambda2=0,tau=0,size_factors=size_factors,norm_y=norm_y,true_clusters=true_clusters,true_disc=true_disc,prefix=pref,dir=dir_name,method=method,disp=disp)  # no penalty
      end=as.numeric(Sys.time())
      list_BIC[aa,2]<-X$BIC
      if(list_BIC[aa,1]==true.K){
        compare_X = X
      }
      print(list_BIC[aa,])
      print(paste("Time:",end-start,"seconds"))
    }
    
  
    sink(file=sprintf("Diagnostics/%s/%s_%s_final%d_order.txt",dir_name,method,disp,ii))
    
    max_k=list_BIC[which.min(list_BIC[,2]),1]
    cat(paste("True order:",true.K,"\n"))
    cat(paste("Optimal order selected:",max_k,"\n"))
    cat("RUN WITH CORRECT ORDER:\n")
    MSE_coefs = sum((compare_X$coefs - sim_coefs[idx,])^2)/(sum(idx)*true.K)
    if(is.null(ncol(init_phi))){
      MSE_phi = sum((init_phi[idx]-compare_X$phi)^2)/(sum(idx)*true.K) # test
    } else{
      MSE_phi = sum(((10/exp(sim_coefs[idx,])^2)-compare_X$phi)^2)/(sum(idx)*true.K) # test
    }
    cat(paste("ARI:",adjustedRandIndex(compare_X$final_clusters,true_clusters),"\n"))
    cat(paste("MSE of true vs discovered coefs:",MSE_coefs,"\n"))
    cat(paste("MSE of true vs discovered phi:",MSE_phi,"\n"))
    cat(paste("% of correctly ID'ed disc genes:",sum(!compare_X$nondiscriminatory==true_disc)/sum(true_disc),"\n"))
    cat(paste("PPs (n x k):\n"))
    write.table(t(compare_X$wts),quote=F,col.names=F)
    
    sink()
    
    pdf(file=sprintf("Diagnostics/%s/%s_%s_final%d_order.pdf",dir_name,method,disp,ii))
    for(c in 1:true.K){
      cl_ids = true_clusters==c
      for(cc in 1:true.K){
        boxplot(compare_X$wts[cc,cl_ids],main=sprintf("Boxplot of PP for subjects of true cl%d being in cl%d",c,cc))
      }
    }
    annotation_col = data.frame(cbind(true_clusters,compare_X$final_clusters))
    colnames(annotation_col)=c("True","Derived")
    rownames(annotation_col)=c(1:ncol(norm_y))
    colnames(norm_y)=c(1:ncol(norm_y))
    annotation_col2 = annotation_col[order(true_clusters),]
    pheatmap(log(norm_y[,order(compare_X$final_clusters)]+0.1),cluster_cols = F,scale="row",annotation_col = annotation_col2)
    dev.off()
    
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
      lambda1_search=1
      lambda2_search=c(0.05,0.1,0.2,0.5,1,1.5,2)
      tau_search=seq(from=0.1,to=0.9,by=0.2)
      
      list_BIC=matrix(0,nrow=length(lambda1_search)*length(lambda2_search)*length(tau_search),ncol=4) # matrix of BIC's: one for each combination of penalty params 
      
      list_BIC[,1]=rep(lambda1_search,each=length(lambda2_search)*length(tau_search))
      list_BIC[,2]=rep(rep(lambda2_search,each=length(tau_search)),times=length(lambda1_search))
      list_BIC[,3]=rep(tau_search,times=length(lambda1_search)*length(lambda2_search))
      
      #search for optimal penalty parameters
      for(aa in 1:nrow(list_BIC)){
        pref = sprintf("grid%d",ii)
        start = as.numeric(Sys.time())
        X<-EM(y=y,k=max_k,tau=list_BIC[aa,3],lambda1=list_BIC[aa,1],lambda2=list_BIC[aa,2],size_factors=size_factors,norm_y=norm_y,true_clusters=true_clusters,true_disc=true_disc,disp=disp,prefix=pref,dir=dir_name,method=method)
        end = as.numeric(Sys.time())
        list_BIC[aa,4]<-X$BIC
        print(list_BIC[aa,])
        print(paste("Time:",end-start,"seconds"))
      }
      
      #store optimal penalty parameters
      max_index<-which(list_BIC[,4]==min(list_BIC[,4]))
      max_tau<-list_BIC[max_index,3]
      max_lambda1<-list_BIC[max_index,1]
      max_lambda2<-list_BIC[max_index,2]
      
      print(paste("Dataset ", ii, "grid search results: ",list_BIC[max_index,]))
      
      if(length(max_index)>1){
        warning("more than one max index")
        max_index<-max_index[1]
        max_tau<-list_BIC[max_index,3]
        max_lambda1<-list_BIC[max_index,1]
        max_lambda2<-list_BIC[max_index,2]
      }
    } else {
      max_tau=0
      max_lambda1=0
      max_lambda2=0
      print(paste("No Penalization"))
    }
    
    # Final run with optimal parameters
    pref = sprintf("final%d",ii)
    start = as.numeric(Sys.time())
    X<-EM(y=y,k=max_k,tau=max_tau,lambda1=max_lambda1,lambda2=max_lambda2,size_factors=size_factors,norm_y=norm_y,true_clusters=true_clusters,true_disc=true_disc,prefix=pref,dir=dir_name,method=method,disp=disp)
    end = as.numeric(Sys.time())
    X$time_elap = end-start
    
    print("Optimal EM run complete")
    
    cls_iClust=NA
    ARI_iClust=NA
    sil_iClust=NA
    # iCluster+
    if(iCluster_compare){
      n.lambda=25
      iCluster_res = sim.iCluster(y[1:200,],true_clusters,ncores=1,n.lambda=n.lambda)
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
      cls_iClust = iCluster_res$cls
      ARI_iClust = iCluster_res$ARI
    } else{iCluster_res=NA}
    
    print("iCluster run complete")
    
    
    # Predictions:
    y_pred=NA
    pred_acc=NA
    true_clusters_pred=NA
    if(max_k == true.K){
      X_pred = X
    
      n_pred = floor(0.1*n)           # simulate data for 10% of original n
      SF_pred = size_factors[1:n_pred]
      
      if(disp=="cluster"){    # check for whether init_phi is of dimension 1
        sim.dat<-simulate_data(n=n_pred,k=true.K,g=g,init_pi=sim_pi,b=sim_coefs,size_factors=SF_pred,distrib=distrib,phi=init_phi) # cluster-wise disp param
      } else{
        sim.dat<-simulate_data_g(n=n_pred,k=true.K,g=g,init_pi=sim_pi,b=sim_coefs,size_factors=SF_pred,distrib=distrib,phi=init_phi) # gene-specific disp param
      }
      y_pred<-sim.dat$y
      z_pred<-sim.dat$z
      true_clusters_pred<-rep(0,times=n_pred)
      
      # same pre-filtering for prediction dataset:
      y_pred = y_pred[idx,]
      
      for(i in 1:n_pred){
        true_clusters_pred[i]<-which(z_pred[,i]==1)
      }
      
      # match cluster ID'sbased on SSE's of true vs estimated coefficients
      SSEs=matrix(0,nrow=true.K,ncol=true.K)
      subs_sim_coefs=sim_coefs[idx,]
      for(c in 1:true.K){
        for(cc in 1:true.K){
          SSEs[c,cc] = sum(abs(subs_sim_coefs[,c]-X_pred$coefs[,cc]))
        }
      }
      new_cl_ids = rep(0,true.K)
      for(c in 1:true.K){
        new_cl_ids[c]=which.min(SSEs[c,])
      }
      true_clusters_pred=new_cl_ids[true_clusters_pred]
      
      
      fit= sim.predict(X_pred,y_pred,SF_pred,true_clusters=true_clusters_pred)
      pred_acc=fit$pred_acc
    }
    
    print("Prediction complete")
    
    pam1 <- function(x,k) list(cluster = pam(t(x),k, cluster.only=TRUE))
    hclusCut <- function(x, k){
      list(cluster = cutree(hclust(as.dist(1-cor(x,method="spearman")),method="average"), k=k))
    }
    
    K_hc=which.max(clusGap(log(norm_y+0.1),FUN=pam1,K.max=7)$Tab[,"gap"])
    K_med=which.max(clusGap(norm_y,FUN=hclusCut,K.max=7)$Tab[,"gap"])
    
    cls_hc = hclusCut(norm_y,K_hc)$cluster
    cls_med = pam1(norm_y,K_med)$cluster
    
    # # Comparisons with Average linkage HC and K-Medoid clustering
    # d<-as.dist(1-cor(norm_y, method="spearman"))  ##Spearman correlation distance w/ log transform##
    # model<-hclust(d,method="average")       # hierarchical clustering
    # 
    # temp_silhc=rep(0,times=6)
    # temp_silmed=rep(0,times=6)
    # temp_clshc=matrix(0,nrow=n,ncol=6)
    # temp_clsmed=matrix(0,nrow=n,ncol=6)
    # for(c in 2:7){
    #   cls_hc <- cutree(model,k=c)
    #   fit = pam(t(norm_y),c)
    #   cls_med = fit$clustering
    #   temp_clshc[,c-1] = cls_hc
    #   temp_clsmed[,c-1] = cls_med
    #   temp_silhc[c-1] = mean(silhouette(cls_hc,d)[,3])
    #   temp_silmed[c-1] = mean(silhouette(cls_med,d)[,3])
    # }
    # hc_id=which.max(temp_silhc)
    # med_id=which.max(temp_silmed)
    # 
    # cls_hc=temp_clshc[,hc_id]
    # cls_med=temp_clsmed[,hc_id]
    # sil_hc=temp_silhc[hc_id]
    # sil_med=temp_silmed[hc_id]
    # K_hc=hc_id+1
    # K_med=med_id+1
    
    cls_EM = X$final_clusters              # in case wrong order is selected. This inputs correct order and outputs final clusters
    
    # d2 = dist(t(norm_y))
    # if(iCluster_compare){
    #   sil_iClust = mean(silhouette(cls_iClust,d2)[,3])
    # }
    # sil_EM = mean(silhouette(cls_EM,d2)[,3])
    
    print(paste("Time:",X$time_elap,"seconds"))
    print(paste("Dataset ",ii,"complete"))
    
    sink()
    
    sil_iClust=NA
    sil_EM=NA
    sil_med=NA
    sil_hc=NA
    
    results=list(X=X,
                 max_k=max_k,
                 max_lambda1=max_lambda1,
                 max_lambda2=max_lambda2,
                 max_tau=max_tau,
                 true_clusters=true_clusters,
                 true_disc=true_disc,
                 iCluster_res=iCluster_res,
                 pred_acc=pred_acc,
                 y_pred=y_pred,
                 true_clusters_pred=true_clusters_pred,
                 cls_hc=cls_hc,
                 cls_med=cls_med,
                 cls_EM=cls_EM,
                 cls_iClust=cls_iClust,
                 ARI_iClust=ARI_iClust,
                 sil_iClust=sil_iClust,
                 sil_EM=sil_EM,
                 sil_med=sil_med,
                 sil_hc=sil_hc,
                 K_hc=K_hc,
                 K_med=K_med)
    return(results)
  }
  
  
  
  temp_ks<-rep(0,times=sim)
  temp_lambda1s<-rep(0,times=sim)
  temp_lambda2s<-rep(0,times=sim)
  temp_taus<-rep(0,times=sim)
  #temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
  #temp_coefs<-list()
  temp_disc<-rep(0,times=sim)
  temp_ARI<-rep(0,times=sim)
  temp_sensitivity<-rep(0,times=sim)
  temp_falsepos<-rep(0,times=sim)
  temp_pred_acc<-rep(0,times=sim)

  ## ADD PARALLELIZATION HERE ##
  cl<-makeCluster(ncores)
  clusterExport(cl=cl,varlist=c(ls(),"EM","EM_run","logsumexpc","soft_thresholding","NB.GOF","simulate_data","simulate_data_g","sim.iCluster","sim.predict","predictions"),envir=environment())
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
  all_iCluster_res <- list(list())
  
  itemp_ARI = rep(0,sim)
  itemp_K = rep(0,sim)
  itemp_lambda = rep(0,sim)
  
  temp_ARI_hc = rep(0,sim)
  temp_ARI_med = rep(0,sim)
  temp_ARI_EM = rep(0,sim)
  temp_ARI_iClust = rep(0,sim)
  
  # temp_sil_hc = rep(0,sim)
  # temp_sil_med = rep(0,sim)
  # temp_sil_EM = rep(0,sim)
  # temp_sil_iClust = rep(0,sim)
  
  temp_K_hc = rep(0,sim)
  temp_K_med = rep(0,sim)
  
  # Summarize results
  for(ii in 1:sim){
    X=par_sim_res[[ii]]$X
    true_clusters=par_sim_res[[ii]]$true_clusters
    true_disc=par_sim_res[[ii]]$true_disc
    
    temp_ks[ii] = par_sim_res[[ii]]$max_k
    temp_lambda1s[ii] = par_sim_res[[ii]]$max_lambda1
    temp_lambda2s[ii] = par_sim_res[[ii]]$max_lambda2
    temp_taus[ii] = par_sim_res[[ii]]$max_tau
    
    temp_pred_acc[ii] <- par_sim_res[[ii]]$pred_acc
    
    temp_disc[ii]<-sum(!X$nondiscriminatory)
    temp_ARI[ii]<-adjustedRandIndex(true_clusters,X$final_clusters)
    
    if(tt>0){
      temp_sensitivity[ii]<-sum(!X$nondiscriminatory[true_disc])/tt                   # tt = # of disc genes simulated (floor(num.disc*g))
    } else {temp_sensitivity[ii]<-NA}
    if(tt<g){
      temp_falsepos[ii]<-sum(!X$nondiscriminatory[!true_disc])/(g-tt)
    } else {temp_falsepos[ii]<-NA}         # take into account genes omitted b/c rowSum of count was <100 (such genes are assumed nondiscriminatory)
    
    all_sim_data[[ii]] = list(par_sim_res[[ii]]$y_pred,par_sim_res[[ii]]$true_clusters_pred)
    if(iCluster_compare){
      all_iCluster_res[[ii]] = par_sim_res[[ii]]$iCluster_res
      itemp_ARI[ii] = all_iCluster_res[[ii]]$ARI
      itemp_K[ii] = all_iCluster_res[[ii]]$K
      itemp_lambda[ii] = all_iCluster_res[[ii]]$lambda
      temp_ARI_iClust[ii] = par_sim_res[[ii]]$ARI_iClust
    }
    
    temp_ARI_hc[ii] = adjustedRandIndex(par_sim_res[[ii]]$cls_hc,true_clusters)
    temp_ARI_med[ii] = adjustedRandIndex(par_sim_res[[ii]]$cls_med,true_clusters)
    temp_ARI_EM[ii] = adjustedRandIndex(par_sim_res[[ii]]$cls_EM,true_clusters)
    
    # temp_sil_hc[ii] = par_sim_res[[ii]]$sil_hc
    # temp_sil_med[ii] = par_sim_res[[ii]]$sil_med
    # temp_sil_EM[ii] = par_sim_res[[ii]]$sil_EM
    # temp_sil_iClust[ii] = par_sim_res[[ii]]$sil_iClust
    
    temp_K_hc[ii] = par_sim_res[[ii]]$K_hc
    temp_K_med[ii] = par_sim_res[[ii]]$K_med
  }
  
  #mean_pi<-colSums(temp_pi)/sim
  #mean_coefs<-Reduce('+',temp_coefs)/sim
  mean_disc<-mean(temp_disc)/g
  mean_ARI<-mean(temp_ARI)
  mean_sensitivity<-mean(temp_sensitivity)
  mean_falsepos<-mean(temp_falsepos)
  
  tab_k = table(temp_ks)
  tab_lambda2 = table(temp_lambda2s)
  tab_tau = table(temp_taus)
  
  #### Final_K, lambda2, tau represent most frequently found K, lambda2, and tau
  # final_k = as.numeric(names(which.max(tab_k)))
  final_k=mean(temp_ks)
  final_lambda2 = as.numeric(names(which.max(tab_lambda2)))
  final_tau = as.numeric(names(which.max(tab_tau)))
  
  #### iCluster+ results
  imean_ARI = mean(itemp_ARI)
  # ifinal_k = as.numeric(names(which.max(table(itemp_K))))
  ifinal_k = mean(itemp_K)
  ifinal_lambda = as.numeric(names(which.max(table(itemp_lambda))))
  
  #### mean prediction accuracy
  mean_pred_acc = mean(temp_pred_acc,na.rm=T)
  
  #### ARI from Average-linkage HC and K-medoid
  mean_ARI_hc = mean(temp_ARI_hc)
  mean_ARI_med = mean(temp_ARI_med)
  mean_ARI_EM = mean(temp_ARI_EM)
  mean_ARI_iClust = mean(temp_ARI_iClust)
  
  # #### Silhouette values comparisons
  # mean_sil_hc = mean(temp_sil_hc)
  # mean_sil_med = mean(temp_sil_med)
  # mean_sil_EM = mean(temp_sil_EM)
  # mean_sil_iClust = mean(temp_sil_iClust)
  mean_sil_hc=NA
  mean_sil_med=NA
  mean_sil_EM=NA
  mean_sil_iClust=NA
  
  # final_K_hc = as.numeric(names(which.max(table(temp_K_hc))))
  # final_K_med = as.numeric(names(which.max(table(temp_K_med))))
  
  final_K_hc = mean(temp_K_hc)
  final_K_med = mean(temp_K_med)
  
  # Store for tabulation:
  
  results<-list(K=final_k,
                all_k=temp_ks,
                lambda2=final_lambda2,
                all_lambda2=temp_lambda2s,
                tau=final_tau,
                all_tau=temp_taus,
                ARI=mean_ARI,
                disc=mean_disc,
                sens=mean_sensitivity,
                falsepos=mean_falsepos,
                all_data=all_data,
                all_iCluster_res=all_iCluster_res,
                imean_ARI=imean_ARI,
                ifinal_k=ifinal_k,
                ifinal_lambda=ifinal_lambda,
                mean_pred_acc=mean_pred_acc,
                all_sim_data=all_sim_data,
                ARI_HC=mean_ARI_hc,
                ARI_med=mean_ARI_med,
                ARI_EM=mean_ARI_EM,
                ARI_iClust=mean_ARI_iClust,
                sil_HC=mean_sil_hc,
                sil_med=mean_sil_med,
                sil_EM=mean_sil_EM,
                sil_iClust=mean_sil_iClust,
                final_K_hc=final_K_hc,
                final_K_med=final_K_med)
  
  return(results)
}
  
  
  
  ###################################
  ###################################
  ###################################
  
  
