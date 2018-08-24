source("sim_EM.R")

sim.EM<-function(true.K, fold.change, num.disc, g, n, 
                 distrib,method="EM",pval_thresh=0.4,filt_method=c("pval","mad"),
                 disp="gene",fixed_parms=T, fixed_coef=6.5,fixed_phi=0.35,
                 nsims=10){
  
  # disp: "gene" or "cluster"
  # low: coef 3.75-3.84, phi 0.13-0.15
  # med: coef 6.59-6.62, phi 0.32-0.38
  # high: coef 7.84-7.85, phi 1.00-1.32
  
  # Fixed phi: scalar for gene-wise, vector of length K for cluster-wise
  
  sim = nsims       # number of sims (set eq to number of cores for now)
  
  dir_name = sprintf("Sim_%d_%d_%d_%f_%f_%s_fixed_%f_%f_%s",n,g,true.K,fold.change,num.disc,distrib,fixed_coef,fixed_phi,filt_method)

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
    
    if(!is.null(ncol(init_phi))){    # check for whether init_phi is of dimension 1
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
    
    
    all_data[[ii]]<-list(y=y,
                         true_clusters=true_clusters,
                         size_factors=size_factors,
                         norm_y=norm_y,
                         true_disc=true_disc
    )
  }
  
  final_Ks = rep(NA,sim)
  
  # Function to run simulation in parallel
  for(ii in 1:sim){
    
    y = all_data[[ii]]$y
    true_clusters = all_data[[ii]]$true_clusters
    norm_y = all_data[[ii]]$norm_y
    true_disc = all_data[[ii]]$true_disc
    
    if(filt_method=="pval"){
      pvals = NB.GOF(y=y,size_factors=size_factors,nsim=1000)
      FDR_pvals = p.adjust(pvals,"fdr")
      # pre-filtering by pval
      filt_ids = (pvals <= pval_thresh)
    } else if(filt_method=="mad"){
      mads = rep(0,g)
      for(j in 1:g){
        mads[j] = mad(log(norm_y[j,]+0.1))
      }
      filt_ids = mads >= quantile(mads,0.75)
    }
    
    y=y[filt_ids,]
    norm_y=norm_y[filt_ids,]
    true_disc=true_disc[filt_ids]
    subs_sim_coefs=sim_coefs[filt_ids,]
    
    if(disp=="gene"){
      subs_init_phi=init_phi[filt_ids]
    } else if(disp=="cluster"){
      subs_init_phi=init_phi[filt_ids,]
    }
    
    # Order selection
    K_search=c(1:7)
    list_BIC=matrix(0,nrow=length(K_search),ncol=2)
    list_BIC[,1]=K_search
    print(paste("Dataset",ii,"Order Selection:"))
    
    for(aa in 1:nrow(list_BIC)){
      pref = sprintf("order%d",ii)
      X<-EM(y=y,k=list_BIC[aa,1],lambda1=0,lambda2=0,tau=0,size_factors=size_factors,norm_y=norm_y,true_clusters=true_clusters,true_disc=true_disc,prefix=pref,dir=dir_name,method=method,disp=disp)  # no penalty
      list_BIC[aa,2]<-X$BIC
      if(list_BIC[aa,1]==true.K){
        compare_X = X
      }
      print(list_BIC[aa,])
      print(paste("Time:",X$time_elap,"seconds"))
    }
    
    sink(file=sprintf("Diagnostics/%s/%s_%s_final%d_order.txt",dir_name,method,disp,ii))
    max_k=list_BIC[which.min(list_BIC[,2]),1]
    cat(paste("True order:",true.K,"\n"))
    cat(paste("Optimal order selected:",max_k,"\n"))
    cat("RUN WITH CORRECT ORDER:\n")
    MSE_coefs = sum((compare_X$coefs - subs_sim_coefs)^2)/(sum(filt_ids)*true.K)
    MSE_phi = sum((subs_init_phi-compare_X$phi)^2)/(sum(filt_ids)*true.K) # test
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
    
    save(compare_X,file=sprintf("Filtering/true_K_%s_%s_run_%d.out",filt_method,disp,ii))
    
    final_Ks[ii] = max_k
  }
  
  return(final_Ks)
}


mad_gene_Ks=sim.EM(true.K=4, fold.change=1, num.disc=0.1, g=1000, n=160, 
       distrib="nb",method="EM",pval_thresh=0.4,filt_method="mad",
       disp="gene",fixed_parms=T, fixed_coef=6.5,fixed_phi=0.35,
       nsims=10)
mad_cl_Ks=sim.EM(true.K=4, fold.change=1, num.disc=0.1, g=1000, n=160, 
       distrib="nb",method="EM",pval_thresh=0.4,filt_method="mad",
       disp="cluster",fixed_parms=T, fixed_coef=6.5,fixed_phi=0.35,
       nsims=10)

pval_gene_Ks=sim.EM(true.K=4, fold.change=1, num.disc=0.1, g=1000, n=160, 
                   distrib="nb",method="EM",pval_thresh=0.4,filt_method="pval",
                   disp="gene",fixed_parms=T, fixed_coef=6.5,fixed_phi=0.35,
                   nsims=10)
pval_cl_Ks=sim.EM(true.K=4, fold.change=1, num.disc=0.1, g=1000, n=160, 
                 distrib="nb",method="EM",pval_thresh=0.4,filt_method="pval",
                 disp="cluster",fixed_parms=T, fixed_coef=6.5,fixed_phi=0.35,
                 nsims=10)
