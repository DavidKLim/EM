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

sourceCpp("M_step.cpp")


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

soft_thresholding=function(alpha,lambda){
  if(abs(alpha)-lambda<0){
    return(0)
  }else{
    return(sign(alpha)*(abs(alpha)-lambda))
  }
}
 
 
# phi.ml <-
#   function(y, mu, n = sum(weights), weights, limit = 10,
#            eps = .Machine$double.eps^0.25,
#            trace = FALSE){
#     lambda = 1E-50            ### change to phi instead of theta? ###
#     score <- function(n, ph, mu, y, w){
#       sum(w*(digamma((1/ph) + y) - digamma(1/ph) + log(1/ph) +
#                1 - log((1/ph) + mu) - (y + (1/ph))/(mu + (1/ph))))*(-1/ph^2) + 2*lambda*ph
#     }
#     info <- function(n, ph, mu, y, w){
#       sum(w*( - trigamma((1/ph) + y) + trigamma((1/ph)) - ph +
#                 2/(mu + (1/ph)) - (y + (1/ph))/(mu + (1/ph))^2))*(1/ph^4) + 2*lambda
#     }
#     if(inherits(y, "lm")) {
#       mu <- y$fitted.values
#       y <- if(is.null(y$y)) mu + residuals(y) else y$y
#     }
#     if(missing(weights)) weights <- rep(1, length(y))
#     #t0 <- n/sum(weights*(y/mu - 1)^2)
#     p0 <- sum(weights*(y/mu - 1)^2)/n
#     it <- 0
#     del <- 1
#     if(trace) message(sprintf("phi.ml: iter %d 'phi = %f'",
#                               it, signif(p0)), domain = NA)
#     while((it <- it + 1) < limit && abs(del) > eps) {
#       p0 <- abs(p0)
#       del <- score(n, p0, mu, y, weights)/(i <- info(n, p0, mu, y, weights))
#       p0 <- p0 + del
#       if(trace) message("phi.ml: iter", it," phi =", signif(p0))
#     }
#     
#     if(p0 < 0) {
#       p0 <- 0
#       warning("estimate truncated at zero")
#       attr(p0, "warn") <- gettext("estimate truncated at zero")
#     }
#     
#     if(it == limit) {
#       warning("iteration limit reached")
#       attr(p0, "warn") <- gettext("iteration limit reached")
#     }
#     attr(p0, "SE") <- sqrt(1/i)
#     res <- list(p0=p0)
#     return(res)
#   }



# M.step<-function(j){
#   
#   
#   beta<-coefs[j,]                                                   # Retrieve beta & theta from
#   theta<-theta_list[[j]]                                            # previous iteration
# 
#   
#   temp<-matrix(0,ncol=(2*k),nrow=maxit_IRLS)    # Temporarily store beta to test for convergence of IRLS
#   dat_j<-dat[dat[,"g"]==j,]                     # subset just the j'th gene
#   
#   
#   for(i in 1:maxit_IRLS){
#     eta <- 
#       if(a==1 & i==1){
#         matrix(rep(beta,times=n),nrow=n,byrow=TRUE)              # first initialization of eta
#       }else if(a>1 & i==1){
#         matrix(rep(beta,times=n),nrow=n,byrow=TRUE) + offset     # Retrieval of eta for IRLS (prev. beta + offset)
#       }else{eta}
#     
#     temp[i,]<-c(beta,phi[j,])
#     
#     for(c in 1:k){
#       
#       dat_jc<-dat_j[dat_j[,"clusts"]==c,]    # subset data for j'th gene, c'th cluster
#       
#       family=negative.binomial(theta=1/phi[j,c])        # can specify family here (plug in updated phi)
#       
#       linkinv<-family$linkinv              # g^(-1) (eta) = mu
#       mu.eta<-family$mu.eta                # g' = d(mu)/d(eta)
#       variance<-family$variance
#       
#       # Estimate beta #
#       mu = linkinv(eta)
#       mu.eta.val = mu.eta(eta)
#       
#       good <- (dat_jc[,"weights"]>0) & (mu.eta.val[,c] != 0)
#       trans_y <- (eta[,c] - offset)[good] + (dat_jc[,"count"][good] - mu[,c][good]) / mu.eta.val[,c][good]    # subtract size factor from transf. y
#       w <- sqrt(dat_jc[,"weights"][good]*mu.eta.val[,c][good]^2/variance(mu[,c])[good])     # weights used in IRLS
#       
#       beta[c] <-
#         if(lambda1 != 0){
#           ((lambda1*((sum(beta)-beta[c]) + (sum(theta[c,])-theta[c,c])))  +  ((1/n)*sum(w*trans_y))) / ((lambda1*(k-1)) + (1/n)*sum(w))
#         } else { beta[c]<-sum(w*trans_y)/sum(w) }
#       
#       #in case beta goes to -inf/+inf: continue with warning
#       if(beta[c]<(-100)){
#         warning(paste("Cluster",c,"Gene",j,"goes to -infinity"))
#         beta[c] = -100
#       }
#       if(beta[c]>100){
#         warning(paste("Cluster",c,"Gene",j,"goes to +infinity"))
#         beta[c] = 100
#       }
#       
#       eta[,c]<-beta[c] + offset      # add back size factors to eta
#       mu[,c]<-linkinv(eta[,c])
#       
#       
#       # Estimate phi #
#   
#       # USING PENALIZED PHI ESTIMATION:::
#       if(all((dat_jc[dat_jc[,"weights"]==1,"count"]-dat_jc[dat_jc[,"weights"]==1,"count"][1])==0)==FALSE){
#         phi[j,c]<- phi.ml(y=dat_jc[,"count"],
#              mu=mu[,c],
#              weights=dat_jc[,"weights"],
#              limit=100,trace=TRUE)$p0
#       } else{phi[j,c]=0}
#       # if condition sets phi = 0 if all observations in cluster are equal
#       ########################################
#       
#     }
#     
#     
#     # update on theta (beta_i - beta_j)
#     for(c in 1:k){
#       for(cc in 1:k){
#         if(abs(theta[c,cc])>=tau){theta[c,cc]<-beta[c]-beta[cc]}        # TLP from Pan paper
#         else{theta[c,cc]<-soft_thresholding(beta[c]-beta[cc],lambda2)}
#       }
#     }
#     
#     # break conditions for IRLS
#     if(i>1){
#       if(sum((temp[i,]-temp[i-1,])^2)<IRLS_tol){   # Sum of Squares of estimated parameters
#         coefs_j<-beta
#         theta_j<-theta
#         temp_j<-temp[1:i,]
#         phi_j<-phi[j,]
#         break
#       }
#     }
#     if(i==maxit_IRLS){
#       coefs_j<-beta
#       theta_j<-theta  
#       temp_j<-temp[1:i,]
#       phi_j<-phi[j,]
#     }
#   }
#   
#   results=list(coefs_j=coefs_j,
#                theta_j=theta_j,
#                temp_j=temp_j,
#                phi_j=phi_j)
#   return(results)
# }




EM_run <- function(y, k,
                   lambda1=0, lambda2=0, tau=0,
                   size_factors=rep(1,times=ncol(y)) ,
                   norm_y=y,
                   true_clusters=NA,
                   init_parms=FALSE,
                   init_coefs=matrix(0,nrow=nrow(y),ncol=k),
                   init_phi=matrix(0,nrow=nrow(y),ncol=k),
                   cls_init,
                   maxit_EM=100){
  
  start_time <- Sys.time()
  
  n<-ncol(y)
  g<-nrow(y)
  
  # adds 0.1 to all y
  y = y+0.1
  
  #no_cores<-detectCores()-1   # for parallel computing
  
  
  # Stopping Criteria
  IRLS_tol = 1E-6
  maxit_IRLS = 50
  EM_tol = 1E-6
  
  # Initialize parameters
  finalwts<-matrix(rep(0,times=k*ncol(y)),nrow=k)
  coefs<-matrix(rep(0,times=g*k),nrow=g)
  pi<-rep(0,times=k)
  phi= matrix(0,nrow=g,ncol=k)    # initial gene-specific dispersion parameters for negative binomial
  # --> Poisson (phi = 0 = 1/theta)
  theta_list <- list()            # temporary to hold all K x K theta matrices across EM iterations
  temp_list <- list()             # store temp to see progression of IRLS
  phi_list <- list()              # store each iteration of phi to see change with each iteration of EM
  
  offset=log(size_factors)
  #offset=rep(0,times=n)            # no offsets
  
  Q<-rep(0,times=maxit_EM)
  
  cls=cls_init

  # initialize weights
  wts<-matrix(rep(0,times=k*ncol(y)),nrow=k)
  for(c in 1:k){
    wts[c,]=(cls==c)^2
  }
  
  # 4.34 gb used
  print(paste("Memory used before EM:",mem_used()))
  
  ########### M / E STEPS #########
  for(a in 1:maxit_EM){
    
    
    if(a==1){         # Initializations for 1st EM iteration
      if(init_parms){
        coefs=init_coefs
        phi=init_phi
      }
      for(j in 1:g){
        if(!init_parms){
          for(c in 1:k){
            coefs[j,c]<-log(mean(as.numeric(y[j,cls==c])))               # Initialize beta
          }
        }
        beta <- coefs[j,]
        theta<-matrix(rep(0,times=k^2),nrow=k)
        for(c in 1:k){
          for(cc in 1:k){
            theta[c,cc]<-soft_thresholding(beta[c]-beta[cc],lambda2)
          }
        }
        theta_list[[j]] <- theta
      }
    }
    
    par_X=rep(list(list()),g)
    
    
    for(j in 1:g){
      # print(paste("j=",j,"a=",a,"k=",k,"lambda1=",lambda1,"lambda2=",lambda2,"tau=",tau,"IRLS_tol=",IRLS_tol,"maxit_IRLS=",maxit_IRLS))
      # print(paste("offset:",offset))
      # print("dat:")
      # print(head(dat))
      # print("y:")
      # print(head(y))
      # print("theta_list[[1]]:")
      # print(theta_list[[1]])
      # print("coefs:")
      # print(head(coefs))
      # print("phi:")
      # print(head(phi))
      if((a>=5 & all(theta_list[[j]]==0))){next}
      sourceCpp("M_step.cpp")
      par_X[[j]] <- M_step(j=j, a=a, y_j=as.integer(y[j,]), all_wts=wts, offset=offset,
                           k=k,theta_list=theta_list,coefs=coefs,phi=phi,
                           lambda1=lambda1,lambda2=lambda2,tau=tau,
                           IRLS_tol=IRLS_tol,maxit_IRLS=maxit_IRLS)
    }
    
    
    for(j in 1:g){
      if((a>=5 & all(theta_list[[j]]==0))){next}
      coefs[j,] <- par_X[[j]]$coefs_j
      theta_list[[j]] <- par_X[[j]]$theta_j
      temp_list[[j]] <- par_X[[j]]$temp_j
      phi[j,] <- par_X[[j]]$phi_j
    }
    
    phi_list[[a]] <- phi
    
    
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
        l[c,i]<-sum(dnbinom(y[,i]-0.1,size=1/phi[,c],mu=exp(coefs[,c] + offset[i]),log=TRUE))    # posterior log like, include size_factor of subj
      }    # subtract out 0.1 that was added earlier
    }
    
    # store and check Q function
    pt1<-(log(pi)%*%rowSums(wts))
    pt2<-sum(wts*l)
    Q[a]<-pt1+pt2
    
    # break condition for EM
    if(a>5){if(abs(Q[a]-Q[a-5])<EM_tol) {
      finalwts<-wts
      break
    }}
    if(a==maxit_EM){finalwts<-wts}
    
    
    # E step
    
    # update on weights
    logdenom = apply(log(pi) + l, 2,logsumexpc)
    for(c in 1:k){
      wts[c,]<-exp(log(pi[c])+l[c,]-logdenom)
    }
    
    # UB and LB on weights
    for(i in 1:n){
      for(c in 1:k){
        if(wts[c,i]<1E-50){
          wts[c,i]=1E-50
        }
        if(wts[c,i]>(1-1E-50)){
          wts[c,i]=1E-50
        }
      }
    }
    
    current_clusters<-rep(0,times=n)
    for(i in 1:n){
      current_clusters[i]<-which.max(wts[,i])
    }
    
    #print(current_clusters)
    print(paste("EM iter",a,"% of cls unchanged (from initial):",sum(current_clusters==cls_init)/n))
    print(paste("Gene1: # of IRLS iterations used in M step:",nrow(temp_list[[1]][rowSums(temp_list[[1]])!=0,])))
    print(paste("coef:",coefs[1,]))
    print(paste("phi:",phi[1,]))
    print(wts)
    
    
  }
  
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
    if(m[j]==1){nondiscriminatory[j]=TRUE}
  }
  pred.nondiscriminatory<-mean(nondiscriminatory)
  
  log_L<-sum(apply(log(pi) + l, 2, logsumexpc))
  BIC=-2*log_L+log(n*g)*(sum(m)+(k-1))         # -2log(L) + log(#obs)*(#parameters estimated). minimum = best. g*k: total params, sum(m): total # of discriminatory genes
  
  end_time <- Sys.time()
  time_elap <- as.numeric(end_time)-as.numeric(start_time)
  
  result<-list(k=k,
               pi=pi,
               coefs=coefs,
               Q=Q[1:a],
               BIC=BIC,
               nondiscriminatory=nondiscriminatory,
               init_clusters=cls,
               final_clusters=final_clusters,
               phi=phi,
               logL=log_L,
               wts=wts,
               time_elap=time_elap,
               lambda1=lambda1,
               lambda2=lambda2,
               tau=tau,
               size_factors=size_factors,
               norm_y=norm_y)
  return(result)
  
}



EM<-function(y, k,
             lambda1=0, lambda2=0, tau=0,
             size_factors=rep(1,times=ncol(y)) ,
             norm_y=y,
             true_clusters=NA,
             init_parms=FALSE,
             init_coefs=matrix(0,nrow=nrow(y),ncol=k),
             init_phi=matrix(0,nrow=nrow(y),ncol=k)){
  
  n = ncol(y)
  g = nrow(y)
  # Initial Clusterings
  
  ## Hierarchical Clustering
  d<-as.dist(1-cor(norm_y, method="spearman"))  ##Spearman correlation distance w/ log transform##
  model<-hclust(d,method="complete")       # hierarchical clustering
  cls_hc <- cutree(model,k=k)
  
  ## K-means Clustering
  cls_km <- kmeans(t(norm_y),k)$cluster

  # Iterate through 2-it EM with each initialization
  all_init_cls <- cbind(cls_hc,cls_km)
  init_cls_BIC <- rep(0,times=ncol(all_init_cls))
  
  for(i in 1:ncol(all_init_cls)){
    ########################## SIMULATION ONLY #############################
    if(is.na(true_clusters)==FALSE){
      all_perms=allPerms(1:k)
      all_clusts=list()
      temp_clust<-rep(0,times=n)
      
      for(iii in 1:nrow(all_perms)){
        for(ii in 1:n){
          temp_clust[ii]<-all_perms[iii,all_init_cls[ii,i]]
        }
        all_clusts[[iii]]<-temp_clust
      }
      
      all_clusts[[nrow(all_perms)+1]]<-all_init_cls[,i]     # contains all permutations of cluster indices
      match_index<-rep(0,times=nrow(all_perms)+1)
      
      for(ii in 1:(nrow(all_perms)+1)){
        match_index[ii]<-mean(true_clusters==all_clusts[[ii]])     # compares each permutation to true --> % of match
      }
      
      all_init_cls[,i]<-all_clusts[[which.max(match_index)]]
    }
    ########################## SIMULATION ONLY #############################
    fit = EM_run(y,k,lambda1,lambda2,tau,size_factors,norm_y,true_clusters=NA,
                              init_parms=init_parms,init_coefs=init_coefs,init_phi=init_phi,
                              cls_init=all_init_cls[,i], maxit_EM=2)
    init_cls_BIC[i] <- fit$BIC
  }
  
  final_init_cls <- all_init_cls[,which.min(init_cls_BIC)]

  if(init_cls_BIC[1]==init_cls_BIC[2]){
    print("Identical initializations")
  } else{ print(colnames(all_init_cls)[which.min(init_cls_BIC)])}
  
  #TESTING RANDOM CLUSTERING
  r_it=100
  rand_inits = matrix(0,nrow=n,ncol=r_it)
  rand_init_BIC = rep(0,r_it)
  for(r in 1:r_it){
    set.seed(r)
    random_cls = sample(1:k,n,replace=TRUE)
    rand_inits[,r] = random_cls
    rand_init_BIC[r] = EM_run(y,k,lambda1,lambda2,tau,size_factors,norm_y,true_clusters=NA,
                              init_parms=init_parms,init_coefs=init_coefs,init_phi=init_phi,
                              cls_init=random_cls, maxit_EM=2)$BIC
    print(rand_init_BIC[r])
  }
  
  rand_final_cls = rand_inits[,which.min(rand_init_BIC)]
  
  # Final run based on optimal initialization
  #results=EM_run(y,k,lambda1,lambda2,tau,size_factors,norm_y,true_clusters,cls_init=final_init_cls,
  #               init_parms=init_parms,init_coefs=init_coefs,init_phi=init_phi)
  results=EM_run(y,k,lambda1,lambda2,tau,size_factors,norm_y,true_clusters,cls_init=rand_final_cls,
                 init_parms=init_parms,init_coefs=init_coefs,init_phi=init_phi)
  return(results)
}

predictions <- function(X,newdata,new_sizefactors){
  # X is the output object from the EM() function
  init_coefs=X$coefs
  init_phi=X$phi
  init_lambda1=X$lambda1
  init_lambda2=X$lambda2
  init_tau=X$tau
  
  ##### EXPORT THIS?? #####
  # library(DESeq2)
  # row_names<-paste("gene",seq(nrow(newdata)))
  # col_names<-paste("subj",seq(ncol(newdata)))
  # cts<-as.matrix(newdata)
  # rownames(cts)<-row_names
  # colnames(cts)<-col_names
  # coldata<-data.frame(matrix(paste("cl"),nrow=ncol(newdata)))
  # rownames(coldata)<-colnames(cts)
  # colnames(coldata)<-"cluster"
  # dds<-DESeqDataSetFromMatrix(countData = cts,
  #                             colData = coldata,
  #                             design = ~ 1)
  # DESeq_dds<-DESeq(dds)
  # init_size_factors<-sizeFactors(DESeq_dds)
  # init_norm_y<-counts(DESeq_dds,normalized=TRUE)
  ##########################
  
  # results = EM(newdata,ncol(init_coefs),init_lambda1,init_lambda2,init_tau,init_size_factors,init_norm_y,
  #              true_clusters=NA,init_parms=TRUE,init_coefs=init_coefs,init_phi=init_phi)
  
  init_size_factors = new_sizefactors
  offset=log(init_size_factors)
  n=ncol(newdata)
  k=ncol(init_coefs)
  
  
  # nb log(f_k(y_i))
  l<-matrix(0,nrow=k,ncol=n)
  for(i in 1:n){
    for(c in 1:k){
      l[c,i]<-sum(dnbinom(newdata[,i],size=1/init_phi[,c],mu=exp(init_coefs[,c] + offset[i]),log=TRUE))    # posterior log like, include size_factor of subj
    }    # subtract out 0.1 that was added earlier
  }
  
  pi=Xtrain$pi
  
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


