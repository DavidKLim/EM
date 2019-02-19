#library(caret)
#library(mclust)
#library(fpc)
#library(amap)
#library(gplots)
#library(parallel)

library(stats)
library(MASS)
library(permute)
library(Rcpp)
library(RcppArmadillo)
library(pryr)
library(mclust)
library(Biobase)

sourceCpp("M_step2.cpp")


logsumexpc=function(v){  
  if(any(is.infinite(v))){ 
    warning("infinite value in v\n") 
  } 
  if(length(v)==1){ return(v[1]) }  
  sv = sort(v, decreasing=TRUE) 
  res = sum(exp(sv[-1] - sv[1]))
  lse = sv[1] + log(1+res)
  lse 
}

soft_thresholding=function(theta,lambda,alpha){
  a=3.7
  if(abs(theta)<=2*lambda*alpha){
    if(abs(theta)<=lambda*alpha){
      return(0)
    } else{
      return(sign(theta)*(abs(theta)-alpha/(1-alpha)))
    }
  }else if(abs(theta)>2*lambda*alpha & abs(theta)<=a*lambda*alpha){
    return(((a-1)*theta - sign(theta)*a*alpha/(1-alpha))/(a-1-1/(lambda*(1-alpha))))
  }else{
    return(theta)
  }
}
 
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


FSCseq<-function(X=NA, y, k,
             lambda=0,alpha=0,
             size_factors=rep(1,times=ncol(y)),
             norm_y=y,
             purity=rep(1,ncol(y)),offsets=rep(0,ncol(y)),              # Offsets: effect of log(covariate) on count #
             true_clusters=NA, true_disc=NA,
             init_parms=FALSE,
             init_coefs=matrix(0,nrow=nrow(y),ncol=k),
             init_phi=matrix(0,nrow=nrow(y),ncol=k),
             init_cls=NA, n_rinits=50, maxit_inits=15,
             disp=c("gene","cluster"),
             method=c("EM","CEM"),
             prefix="", dir="NA"){
  
  # y: raw counts
  # k: #clusters
  # size_factors: SF's derived from DESeq2
  # norm_y: counts normalized for sequencing depth by DESeq2
  # purity: input custom values to weight E step weights by purity (DISABLED)
  # offsets: option to include additional offsets. must be length n 
  # true_clusters: if applicable. For diagnostics tracking of ARI
  # true_disc: if applicable. For diagnostics tracking of disc/nondisc genes
  # init_parms: TRUE if initial coefficient estimates/dispersion estimates are input
  # init_coefs & init_phi: Initial estimates, if applicable
  # disp = c(gene, cluster), depending on whether dispersions are gene-specific or cluster-specific
  # init_cls: Initial clustering
  # n_rinits: Number of initial clusterings searched with maxit=15. More initializations = more chance to attain global max
  
  if(is.na(X)){
    warning("No covariates specified. Running cluster-specific intercept-only model.")
    } else{
  if (class(X) != "matrix") {
    tmp <- try(X <- model.matrix(~0+., data=X), silent=TRUE)
    if (class(tmp)[1] == "try-error") stop("X must be a matrix or able to be coerced to a matrix")
  }
  if (storage.mode(X)=="integer") storage.mode(X) <- "double"
  if (ncol(y) != nrow(X)) stop("X and y do not have the same number of observations")
  if (any(is.na(y)) | any(is.na(X))) stop("Missing data (NA's) detected.  Take actions (e.g., removing cases, removing features, imputation) to eliminate missing data before passing X and y")
    }
  
  if(alpha==1){
    stop("alpha must be less than 1; choose a smaller number instead")
  } else if(alpha < 0 | alpha > 1){
    stop("alpha not within range [0,1)")
  }
  if(lambda<0){
    stop("lambda must be greater than 0")
  }
  
  if(method=="EM"){
    CEM=F
  } else if(method=="CEM"){
    CEM=T
  } else{
    stop("method must be 'EM' or 'CEM'.")
  }
  
  diag_file = sprintf("Diagnostics/%s/%s_%s_%s_%d_%f_%f.txt",dir,method,disp,prefix,k,lambda,alpha)
  dir.create(sprintf("Diagnostics/%s",dir))
  sink(file=diag_file)
  
  n = ncol(y)
  g = nrow(y)
  
  cat(paste(sprintf("n=%d, g=%d, k=%d, l=%f, alph=%f, ",n,g,k,lambda,alpha),"\n"))
  cat("True clusters:\n")
  write.table(true_clusters,quote=F,col.names=F)
  
  init_Tau=1
  if(CEM){
    init_Tau=sqrt(g)
  }
  if(k==1){
    init_cls=rep(1,n)
  }
  if(all(is.na(init_cls)) & k>1){
    # Initial Clusterings
    ## Hierarchical Clustering
    d<-as.dist(1-cor(norm_y, method="spearman"))  ##Spearman correlation distance w/ log transform##
    model<-hclust(d,method="average")       # hierarchical clustering
    cls_hc <- cutree(model,k=k)
    
    ## K-means Clustering
    cls_km <- kmeans(t(log(norm_y+0.1)),k)$cluster
    
    #TESTING RANDOM CLUSTERING
    
    r_it=n_rinits
    rand_inits = matrix(0,nrow=n,ncol=r_it)
    
    for(r in 1:r_it){
      set.seed(r)
      rand_inits[,r] = sample(1:k,n,replace=TRUE)
      while(sum(1:k %in% rand_inits[,r]) < k){
        rand_inits[,r] = sample(1:k,n,replace=TRUE)
      }
    }
    colnames(rand_inits) = paste("rand",c(1:r_it),sep="")
    
    # Iterate through 2-it EM with each initialization
    all_init_cls <- cbind(cls_hc,cls_km,rand_inits)
    init_cls_BIC <- rep(0,times=ncol(all_init_cls))
    
    all_fits = list()
    
    for(i in 1:ncol(all_init_cls)){
      if(!is.na(X)){
        if(abs(cor(all_init_cls[,i],X))==1){
        init_cls_BIC[i] = .Machine$double.xmax
        next
        }
      }
      
      cat(paste("INITIAL CLUSTERING:",colnames(all_init_cls)[i],"\n"))
      
      fit = EM_run(X,y,k,lambda,alpha,size_factors,norm_y,purity,offsets,true_clusters,true_disc,
                   init_parms=init_parms,init_coefs=init_coefs,init_phi=init_phi,disp=disp,
                   cls_init=all_init_cls[,i], CEM=CEM,init_Tau=init_Tau,maxit_EM=maxit_inits)
      all_fits [[i]] = fit
      init_cls_BIC[i] <- fit$BIC
    }
    
    fit_id = which.min(init_cls_BIC)
    
    cat("FINAL INITIALIZATION:\n")
    cat(paste(colnames(all_init_cls)[fit_id],"\n"))
    init_cls = all_init_cls[,fit_id[1]]
    if(CEM){
      cat("FINAL TAU VALUE:\n")
      cat(paste(init_Tau,"\n"))
    }
  }
  
  results=EM_run(X,y,k,lambda,alpha,size_factors,norm_y,purity,offsets,true_clusters,true_disc,
                 init_parms=init_parms,init_coefs=init_coefs,init_phi=init_phi,disp=disp,
                 cls_init=init_cls, CEM=CEM,init_Tau=init_Tau,maxit_EM=100)
  
  sink()
  return(results)
}





EM_run <- function(X=NA, y, k,
                   lambda=0,alpha=0,
                   size_factors=rep(1,times=ncol(y)) ,
                   norm_y=y,
                   purity=rep(1,ncol(y)),offsets=rep(0,ncol(y)),
                   true_clusters=NA, true_disc=NA,
                   init_parms=FALSE,
                   init_coefs=matrix(0,nrow=nrow(y),ncol=k),
                   init_phi=matrix(0,nrow=nrow(y),ncol=k),
                   disp,cls_init,
                   CEM=F,init_Tau=1,
                   maxit_EM=100){
  
  start_time <- Sys.time()
  
  n<-ncol(y)         # number of samples
  g<-nrow(y)         # number of genes
  p<-ncol(X)         # number of covariates
  if(is.null(p)){p=0}
  
  cl_X = matrix(0,nrow=k*n,ncol=k)
  ident_k = diag(k)
  for(i in 1:k){
    cl_X[((i-1)*n+1):(i*n),] = matrix(rep(ident_k[i,],n),ncol=k,nrow=n,byrow=T)
  }
  
  if(sum(is.na(X))>0){
    XX=cl_X
    covars=F
  } else{
    XX = do.call("rbind", replicate(k, X,simplify=FALSE))
    XX = cbind(cl_X,XX)
    covars=T
  }
  
  # adds 0.1 to all y
  y = y+0.1
  
  #no_cores<-detectCores()-1   # for parallel computing
  
  if(disp=="gene"){
    cl_phi=0
  } else if(disp=="cluster"){
    cl_phi=1
  }
  
  # Stopping Criteria
  IRLS_tol = 1E-12     # for both phi and beta
  maxit_IRLS = 50
  EM_tol = 1E-12
  
  # Initialize parameters
  finalwts<-matrix(rep(0,times=k*ncol(y)),nrow=k)
  coefs<-matrix(rep(0,times=(g*(k+p))),nrow=g)
  pi<-rep(0,times=k)
  phi= matrix(0,nrow=g,ncol=k)    # initial gene-specific dispersion parameters for negative binomial
  # --> Poisson (phi = 0 = 1/theta)
  theta_list <- list()            # temporary to hold all K x K theta matrices across EM iterations
  temp_list <- list()             # store temp to see progression of IRLS
  phi_list <- list()              # store each iteration of phi to see change with each iteration of EM
  coefs_list <- list()
  LFCs = matrix(0,nrow=g,ncol=maxit_EM)
  
  offset=log2(size_factors) + offsets
  #offset=rep(0,times=n)            # no offsets
  
  Q<-rep(0,times=maxit_EM)
  
  cls=cls_init

  # Initialize weights
  wts<-matrix(rep(0,times=k*ncol(y)),nrow=k)
  for(c in 1:k){
    wts[c,]=(cls==c)^2
  }
  
  # For use in CEM in E step #
  Tau = init_Tau
  if(CEM){cat(paste("Tau =",Tau,"\n"))}
  
  phi_g = rep(0,times=g)
  DNC=0
  
  disc_ids_list = list()
  disc_ids=rep(T,g)
  
  diff_phi=matrix(0,nrow=maxit_EM,ncol=g)
  est_phi=rep(1,g)                          # 1 for true, 0 for false
  
  ########### M / E STEPS #########
  for(a in 1:maxit_EM){
    EMstart= as.numeric(Sys.time())
    
    if(a==1){         # Initializations for 1st EM iteration
      if(init_parms){
        coefs=init_coefs
        if(cl_phi==0){
          phi_g=init_phi
        } else{
          phi=init_phi
        }
      }
      for(j in 1:g){
        if(!init_parms){
          fit = glm(as.integer(rep(y[j,],k))~0+XX+offset(rep(offset,k)),family=poisson(),weights=c(t(wts)))
          coefs[j,] = fit$coefficients
          if(covars){
            for(c in (k+1):(k+p)){
              if(is.na(coefs[j,c])){coefs[j,c]=0}
            }
          }
        }
        beta <- coefs[j,]
        theta<-matrix(rep(0,times=k^2),nrow=k)
        for(c in 1:k){
          for(cc in 1:k){
            theta[c,cc]<-soft_thresholding(beta[c]-beta[cc],lambda,alpha)
          }
        }
        theta_list[[j]] <- theta
      }
    }
    
    par_X=rep(list(list()),g)
    
    #sourceCpp("M_step2.cpp")         # TESTING THE NEW M_step FUNCTION
    Mstart=as.numeric(Sys.time())
    for(j in 1:g){
      if(Tau<=1 & a>6){if(Reduce("+",disc_ids_list[(a-6):(a-1)])[j]==0){next}}
      y_j = as.integer(y[j,])
      par_X[[j]] <- M_step(X=XX, p=p, j=j, a=a, y_j=y_j, all_wts=wts, offset=rep(offset,k),
                           k=k,theta=theta_list[[j]],coefs_j=coefs[j,],phi_j=phi[j,],cl_phi=cl_phi,phi_g=phi_g[j],est_phi=est_phi[j],
                           lambda=lambda,alpha=alpha,
                           IRLS_tol=IRLS_tol,maxit_IRLS=maxit_IRLS #,fixed_phi = phis
                           )
    }
    Mend=as.numeric(Sys.time())
    
    
    for(j in 1:g){
      if(Tau<=1 & a>6){if(Reduce("+",disc_ids_list[(a-6):(a-1)])[j]==0){next}}
      coefs[j,] <- par_X[[j]]$coefs_j
      theta_list[[j]] <- par_X[[j]]$theta_j
      temp_list[[j]] <- par_X[[j]]$temp_j
      if(cl_phi==1){
        phi[j,] <- par_X[[j]]$phi_j
      } else if(cl_phi==0){
        phi_g[j] <- (par_X[[j]]$phi_j)[1]
        phi[j,] = rep(phi_g[j],k)
      }
      
      # Ad hoc averaging of cluster log means when nondiscriminatory
      ids = list()
      n_k = rowSums(wts)
      for(c in 1:k){
        ids[[c]]=which(theta_list[[j]][c,]==0)
        coefs[j,ids[[c]]] = rep(sum(n_k[ids[[c]]]*coefs[j,ids[[c]]])/sum(n_k[ids[[c]]]),times=length(ids[[c]]))               # weighted (by # in each cl) average
      }
    }
    # Marker of all nondisc genes (T for disc, F for nondisc)
    for(j in 1:g){
      disc_ids[j]=any(theta_list[[j]]!=0)
    }
    cat(paste("Disc genes:",sum(disc_ids),"of",g,"genes.\n"))
    disc_ids_list[[a]] = disc_ids
    
    if(cl_phi==1){
      phi_list[[a]] <- phi
    } else if(cl_phi==0){
      phi_list[[a]] <- phi_g
    }
    
    if(a>6){
      for(j in 1:g){
        if(cl_phi==1){
          diff_phi[a,j]=mean(abs(phi_list[[a]][j,]-phi_list[[a-5]][j,])/phi_list[[a-5]][j,])
        } else if(cl_phi==0){
          diff_phi[a,j]=abs(phi_list[[a]][j]-phi_list[[a-5]][j])/phi_list[[a-5]][j]
        }
        
        if(is.infinite(diff_phi[a,j]) | is.na(diff_phi[a,j])){
          diff_phi[a,j]=1
        }
        
        if(diff_phi[a,j]<0.01){
          est_phi[j]=0
        } else{
          est_phi[j]=1
        }
      }
      cat(paste("Avg % diff in phi est (across 5 its) gene 1 = ",diff_phi[a,1],"\n"))
    }
    
    coefs_list[[a]] = coefs
    LFCs[,a] = (rowMax(coefs)-rowMin(coefs))/(k-1)
    
    
    
    # update on pi_hat, and UB & LB on pi
    for(c in 1:k){
      pi[c]=mean(wts[c,])
      if(pi[c]<1E-6){
        warning(paste("cluster proportion", c, "close to 0"))
        pi[c]=1E-6
      } # lowerbound for pi
      if(pi[c]>(1-1E-6)){
        warning(paste("cluster proportion", c, "close to 1"))
        pi[c]=(1-1E-6)
      } # upperbound for pi
    }
    
    # nb log(f_k(y_i))
    l<-matrix(rep(0,times=k*n),nrow=k)
    for(i in 1:n){
      for(c in 1:k){
        if(covars){
          covar_coefs = matrix(coefs[,-(1:k)],ncol=p)
          cov_eff = apply(covar_coefs,1,sum)
        } else {cov_eff=rep(0,g)}
        
        if(cl_phi==1){
          l[c,i]<-sum(dnbinom(y[,i]-0.1,size=1/phi[,c],mu=2^(coefs[,c] + cov_eff + offset[i]),log=TRUE))    # posterior log like, include size_factor of subj
        } else if(cl_phi==0){
          l[c,i]<-sum(dnbinom(y[,i]-0.1,size=1/phi_g,mu=2^(coefs[,c] + cov_eff + offset[i]),log=TRUE))
        }
      }    # subtract out 0.1 that was added earlier
    }
    
    # store and check Q function
    pt1<-(log(pi)%*%rowSums(wts))
    pt2<-sum(wts*l)
    Q[a]<-pt1+pt2
    
    # break condition for EM
    if(Tau==1 & a>5){if(abs((Q[a]-Q[a-5])/Q[a])<EM_tol) {
      finalwts<-wts
      break
    }}
    if(a==maxit_EM){
      finalwts<-wts
      warning("Reached max iterations.")
      DNC=1
      }
  
    
    # E step
    
    prev_clusters<-rep(0,times=n)
    for(i in 1:n){
      prev_clusters[i]<-which.max(wts[,i])
    }
    if(a==1 & !is.null(true_clusters) & !is.na(true_clusters)){
      cat(paste("Initial ARI:",adjustedRandIndex(prev_clusters,true_clusters),"\n"))
    }
    
    if(!CEM){
      # update on weights
      logdenom = apply(log(pi) + l, 2,logsumexpc)
      for(c in 1:k){
        wts[c,]<-exp(log(pi[c])+l[c,]-logdenom)
      }
    } else if(CEM){
      # CEM update on weights
      logdenom = apply((1/Tau)*(log(pi)+l),2,logsumexpc)
      for(c in 1:k){
        wts[c,]<-exp((1/Tau)*(log(pi[c])+l[c,])-logdenom)
      }
      if(Tau>1){
        Tau = 0.9*Tau
      } else{ Tau=1 }     # after Tau hits 1 --> fix at 1
    }
    
    # UB and LB on weights
    for(i in 1:n){
      for(c in 1:k){
        if(is.na(wts[c,i])){
          wts[c,i]=1E-50
        } else if(wts[c,i]<1E-50){
          wts[c,i]=1E-50
        } else if(wts[c,i]>(1-1E-50)){
          wts[c,i]=1-1E-50
        }
      }
    }

    if(CEM){
      # CEM
      draw_wts=wts                 # initialize
      for(i in 1:n){
        set.seed(i)
        draw_wts[,i] = rmultinom(1,1,wts[,i])
      }
      seed_mult=1
      while(any(rowSums(draw_wts)==0)){
        cat("Drawing again",seed_mult,"\n")
        for(i in 1:n){
          set.seed(seed_mult*n+i)
          for(c in 1:k){
            if(wts[c,i]<=(1E-50*10^seed_mult) & seed_mult<=48){
              wts[c,i]=1E-50*10^seed_mult
            } else if(wts[c,i]>=(1-(1E-50*10^seed_mult)) & seed_mult<=48){
              wts[c,i]=1-1E-50*10^seed_mult
            }
          }
          draw_wts[,i] = rmultinom(1,1,wts[,i])
        }
        seed_mult=seed_mult+1
        if(seed_mult>250){
          draw_wts[,n]=rep(1/k,k)
          break
          }
      }
      wts=draw_wts
    }          # Keep drawing until at least one in each cluster
    
    # # Weighting by purity : found that this makes no difference
    # for(i in 1:n){
    #   wts[,i]=wts[,i]*purity[i]
    # }      # default: purity = 1 --> wts stay the same.
    
    # Diagnostics Tracking
    current_clusters<-rep(0,times=n)
    for(i in 1:n){
      current_clusters[i]<-which.max(wts[,i])
    }
    
    #print(current_clusters)
    cat(paste("EM iter",a,"% of cls unchanged (from previous):",sum(current_clusters==prev_clusters)/n,"\n"))
    if(!is.null(true_clusters) & !is.na(true_clusters)){cat(paste("ARI =",adjustedRandIndex(true_clusters,current_clusters),"\n"))}
    cat(paste("Cluster proportions:",pi,"\n"))
    if(sum(is.na(true_disc))==0){
      if(sum(true_disc)==0){
        disc_gene=1
        cat("No discriminatory genes. Printing Gene1 instead\n")
        } else{ disc_gene = which(true_disc^2==1)[1] }
      if(sum(true_disc)==length(true_disc)){
        nondisc_gene=2
        cat("No nondiscriminatory genes. Printing Gene1 instead\n")
      } else{ nondisc_gene = which(true_disc^2==0)[1] }
      cat(paste("Disc Gene",disc_gene,": # of IRLS iterations used in M step:",nrow(temp_list[[disc_gene]][rowSums(temp_list[[disc_gene]])!=0,]),"\n"))
      cat(paste("coef:",coefs[disc_gene,],"\n"))
      if(cl_phi==1){
        cat(paste("phi:",phi[disc_gene,],"\n"))
      } else if(cl_phi==0){
        cat(paste("phi:",phi_g[disc_gene],"\n"))
      }
      cat(paste("Nondisc Gene",nondisc_gene,": # of IRLS iterations used in M step:",nrow(temp_list[[nondisc_gene]][rowSums(temp_list[[nondisc_gene]])!=0,]),"\n"))
      cat(paste("coef:",coefs[nondisc_gene,],"\n"))
      if(cl_phi==1){
        cat(paste("phi:",phi[nondisc_gene,],"\n"))
      } else if(cl_phi==0){
        cat(paste("phi:",phi_g[nondisc_gene],"\n"))
      }
    } else{
      cat(paste("Gene1: # of IRLS iterations used in M step:",nrow(temp_list[[1]][rowSums(temp_list[[1]])!=0,]),"\n"))
      cat(paste("coef:",coefs[1,],"\n"))
      if(cl_phi==1){
        cat(paste("phi:",phi[1,],"\n"))
      } else if(cl_phi==0){
        cat(paste("phi:",phi_g[1],"\n"))
      }
    }
    cat(paste("Samp1: PP:",wts[,1],"\n"))
    EMend = as.numeric(Sys.time())
    cat(paste("Time elapsed:",EMend-EMstart,"seconds\n"))
    cat("-------------------------------------\n")
    
  }
  
  cat("-------------------------------------\n")
  num_warns=length(warnings())
  
  final_clusters<-rep(0,times=n)
  for(i in 1:n){
    final_clusters[i]<-which.max(finalwts[,i])
  }
  
  m<-rep(0,times=g)
  nondiscriminatory=rep(FALSE,times=g)
  for(j in 1:g){
    m_row=rep(0,k)
    for(c in 1:k){
      m_row[c] <- sum(theta_list[[j]][c,]!=0) + 1         # of parameters estimated
    }
    m[j]=min(m_row)
    if(sum(coefs[j,]==-30)>1){
      m[j]=m[j]+(sum(coefs[j,]==-30)-1)
    }
    if(m[j]==1){nondiscriminatory[j]=TRUE}
  }
  if(lambda*alpha==0){
    m=rep(k,g)
  }
  
  num_est_coefs = sum(m)
  num_est_params = 
    if(cl_phi==1){
      2*sum(m)+(k-1)+p*g
      } else{ sum(m)+(k-1)+g+p*g }
  
  log_L<-sum(apply(log(pi) + l, 2, logsumexpc))
  BIC = -2*log_L + log(n)*num_est_params
  
  cat(paste("total # coefs estimated =",num_est_coefs,"\n"))
  cat(paste("total # params estimated =",num_est_params,"\n"))
  cat(paste("-2log(L) =",-2*log_L,"\n"))
  cat(paste("log(n) =",log(n),"\n"))
  cat(paste("BIC =",BIC,"\n"))
  
  disc_stats=cbind(m,(!nondiscriminatory)^2,disc_ids^2)
  colnames(disc_stats) = c("#params","disc","disc_ids")
  write.table(head(disc_stats,n=10),quote=F)
  write.table(tail(disc_stats,n=10),quote=F)
  cat("-------------------------------------\n")
  cat("Coefs:\n")
  write.table(head(coefs,n=10),quote=F)
  write.table(tail(coefs,n=10),quote=F)
  cat("-------------------------------------\n")
  cat("Phi:\n")
  if(cl_phi==1){
    write.table(head(phi,n=10),quote=F)
    write.table(tail(phi,n=10),quote=F)
  } else if(cl_phi==0){
    write.table(head(phi_g,n=10),quote=F)
    write.table(tail(phi_g,n=10),quote=F)
  }
  cat("-------------------------------------\n")
  
  
  end_time <- Sys.time()
  time_elap <- as.numeric(end_time)-as.numeric(start_time)
  
  #cat("m (number of params est'ed per gene):\n")
  #write.table(t(m),quote=F,col.names=F,row.names=F)
  #cat("BIC:",BIC,"\n")
  if(cl_phi==0){
    phi = phi_g
  }
  
  result<-list(k=k,
               pi=pi,
               coefs=coefs,
               Q=Q[1:a],
               BIC=BIC,
               nondiscriminatory=nondiscriminatory,
               init_clusters=cls_init,
               final_clusters=final_clusters,
               phi=phi,
               logL=log_L,
               wts=wts,
               time_elap=time_elap,
               lambda=lambda,
               alpha=alpha,
               size_factors=size_factors,
               norm_y=norm_y,DNC=DNC,LFCs=LFCs,disc_ids_list=disc_ids_list)
  return(result)
  
}

predictions <- function(X,newdata,new_sizefactors,purity=rep(1,ncol(newdata)),offsets=rep(0,ncol(newdata))){
  # X: Output of EM
  # newdata: Data to perform prediction on
  # new_sizefactors: SF's of new data
  # purity: Custom purity values can be input and adjusted for (NOT AVAILABLE YET)
  # offsets: Additional offsets per sample can be incorporated
  
  # X is the output object from the EM() function
  init_coefs=X$coefs
  init_phi=X$phi
  init_lambda=X$lambda
  init_alpha=X$alpha
  
  
  cl_phi=!is.null(dim(init_phi))  # dimension of phi is null when gene-wise (vector)
  
  init_size_factors = new_sizefactors
  offset=log2(init_size_factors) + offsets
  n=ncol(newdata)
  k=ncol(init_coefs)
  
  
  # nb log(f_k(y_i))
  l<-matrix(0,nrow=k,ncol=n)
  for(i in 1:n){
    for(c in 1:k){
      if(cl_phi){
        l[c,i]<-sum(dnbinom(newdata[,i],size=1/init_phi[,c],mu=2^(init_coefs[,c] + offset[i]),log=TRUE))    # posterior log like, include size_factor of subj
      } else if(!cl_phi){
        l[c,i]<-sum(dnbinom(newdata[,i],size=1/init_phi,mu=2^(init_coefs[,c] + offset[i]),log=TRUE))
      }
    }    # subtract out 0.1 that was added earlier
  }
  
  pi=X$pi
  
  # E step
  # Estimate weights
  wts = matrix(0,nrow=k,ncol=n)
  logdenom = apply(log(pi) + l, 2,logsumexpc)
  for(c in 1:k){
    wts[c,]<-exp(log(pi[c])+l[c,]-logdenom)
  }
  
  final_clusters<-rep(0,times=n)
  for(i in 1:n){
    final_clusters[i]<-which.max(wts[,i])
  }
  
  return(list(final_clusters=final_clusters,wts=wts))
}
  