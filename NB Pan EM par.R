library(caret)
library(mclust)
library(stats)
library(MASS)
library(fpc)
library(permute)
library(amap)
library(gplots)
library(parallel)



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

theta.ml2 <-
  function(y, mu, n = sum(weights), weights, limit = 10,
           eps = .Machine$double.eps^0.25,
           trace = FALSE,use_ml=1){
    score <- function(n, th, mu, y, w)
      sum(w*(digamma(th + y) - digamma(th) + log(th) +
               1 - log(th + mu) - (y + th)/(mu + th)))
    info <- function(n, th, mu, y, w)
      sum(w*( - trigamma(th + y) + trigamma(th) - 1/th +
                2/(mu + th) - (y + th)/(mu + th)^2))
    if(inherits(y, "lm")) {
      mu <- y$fitted.values
      y <- if(is.null(y$y)) mu + residuals(y) else y$y
    }
    if(missing(weights)) weights <- rep(1, length(y))
    t0 <- n/sum(weights*(y/mu - 1)^2)
    it <- 0
    del <- 1
    if(trace) message(sprintf("theta.ml: iter %d 'theta = %f'",
                              it, signif(t0)), domain = NA)
    while((it <- it + 1) < limit && abs(del) > eps && use_ml==1) {
      t0 <- abs(t0)
      del <- score(n, t0, mu, y, weights)/(i <- info(n, t0, mu, y, weights))
      if(del < (-t0) || is.na(t0) || is.na(del)){
        warning("Theta goes from (+) to (-). Last iteration", it," theta =", signif(t0),". Using method of moments instead")
        use_ml <- 0
        break                 # if the delta is changing the sign of t0 from + to -, then break (keep last iteration of t0)
      }
      t0 <- t0 + del
      if(trace) message("theta.ml: iter", it," theta =", signif(t0))
    }
    
    if(use_ml==0){t0 = theta.mm(y=y,mu=mu,dfr=n-1,weights=weights)}
    
    if(t0 < 0) {
      t0 <- 0
      warning("estimate truncated at zero")
      attr(t0, "warn") <- gettext("estimate truncated at zero")
    }
    if(it == limit) {
      warning("iteration limit reached")
      attr(t0, "warn") <- gettext("iteration limit reached")
    }
    attr(t0, "SE") <- sqrt(1/i)
    res <- list(t0=t0,
                use_ml=use_ml)
    return(res)
  }

M.step<-function(j){
  beta<-rep(0,times=k)
  
  if(a==1){
    for(c in 1:k){
      beta[c]<-log(mean(as.numeric(y[j,cls==c])))               # Initialize beta
    }
    
    theta<-matrix(rep(0,times=k^2),nrow=k)
    for(c in 1:k){
      for(cc in 1:k){
        if(theta[c,cc]>=tau){theta[c,cc]<-beta[c]-beta[cc]}             # Initialize theta
        else{theta[c,cc]<-soft_thresholding(beta[c]-beta[cc],lambda2)   # 
        }
      }
    }
  } else {
    beta<-coefs[j,]                                                   # Retrieve beta & theta from
    theta<-theta_list[[j]]                                            # previous iteration
  }
  
  
  temp<-matrix(0,ncol=(2*k),nrow=maxit_IRLS)    # Temporarily store beta to test for convergence of IRLS
  dat_j<-dat[dat[,"g"]==j,]                                  # subset just the j'th gene
  
  
  for(i in 1:maxit_IRLS){
    eta <- 
      if(a==1 & i==1){
        matrix(rep(beta,times=n),nrow=n,byrow=TRUE)               # first initialization of eta
      }else if(a>1 & i==1){
        matrix(rep(beta,times=n),nrow=n,byrow=TRUE) + offset     # Retrieval of eta for IRLS (prev. beta + offset)
      }else{eta}
    
    
    temp[i,]<-c(beta,phi[j,])
    
    for(c in 1:k){
      
      dat_jc<-dat_j[dat_j[,"clusts"]==c,]    # subset j'th gene, c'th cluster
      
      #phi[j,c]<-glm.nb((dat_jc[,"count"]-0.1) ~ 1 + offset(offset), weights=dat_jc[,"weights"])$theta   #glm.nb() theta
      
      family=negative.binomial(theta=1/phi[j,c])        # can specify family here (plug in updated phi)
      
      linkinv<-family$linkinv              # g^(-1) (eta) = mu
      mu.eta<-family$mu.eta                # g' = d(mu)/d(eta)
      variance<-family$variance
      
      mu = linkinv(eta)
      mu.eta.val = mu.eta(eta)
      
      
      good <- (dat_jc[,"weights"]>0) & (mu.eta.val[,c] != 0)
      
      trans_y <- (eta[,c] - offset)[good] + (dat_jc[,"count"][good] - mu[,c][good]) / mu.eta.val[,c][good]    # subtract size factor from transf. y
      
      w <- sqrt(dat_jc[,"weights"][good]*mu.eta.val[,c][good]^2/variance(mu[,c])[good])     # weights used in IRLS
      
      beta[c] <-
        if(lambda1 != 0){
          ((lambda1*((sum(beta)-beta[c]) + (sum(theta[c,])-theta[c,c])))  +  ((1/n)*sum(w*trans_y))) / ((lambda1*(k-1)) + (1/n)*sum(w))
        } else { beta[c]<-sum(w*trans_y)/sum(w) }
      
      
      if(beta[c]<(-100)){
        warning(paste("Cluster",c,"Gene",j,"goes to -infinity"))
        beta[c] = -100
      }
      if(beta[c]>100){
        warning(paste("Cluster",c,"Gene",j,"goes to +infinity"))
        beta[c] = 100
      }
      
      eta[,c]<-beta[c] + offset      # add back size factors to eta
      mu[,c]<-linkinv(eta[,c])
      
      
      # Calculate phi = 1/theta #
      
      #### Maximum Likelihood Estimation ####
      if(all((dat_jc[dat_jc[,"weights"]==1,"count"]-dat_jc[dat_jc[,"weights"]==1,"count"][1])==0)==FALSE){
          fit <- theta.ml2(y = dat_jc[,"count"]-0.1,
                           mu = mu[,c],
                           weights = dat_jc[,"weights"],
                           limit=10,trace=FALSE,use_ml=phi_use_ml[j,c])
          phi[j,c] <- 1/fit$t0    # this bypasses error when all counts in cluster are identical or
                                # there is just one subject in cluster (this would be 0 disp anyway)
          phi_use_ml[j,c] = fit$use_ml
        } else {phi[j,c]=0}
      ########################################
      
      
    }
    
    
    # update on theta (beta_i - beta_j)
    for(c in 1:k){
      for(cc in 1:k){
        if(abs(theta[c,cc])>=tau){theta[c,cc]<-beta[c]-beta[cc]}                      # gTLP from Pan paper
        else{theta[c,cc]<-soft_thresholding(beta[c]-beta[cc],lambda2)}           # 
      }
    }
    # break conditions for IRLS
    if(i>1){
      if(sum((temp[i,]-temp[i-1,])^2)<IRLS_tol){  # convergence of just coefs?? or phi too?
        coefs_j<-beta                                 # reached convergence
        theta_j<-theta
        temp_j<-temp[1:i,]
        phi_j<-phi[j,]
        phi_use_ml_j<-phi_use_ml[j,]
        break
      }
    }
    if(i==maxit_IRLS){
      coefs_j<-beta
      theta_j<-theta  
      temp_j<-temp[1:i,]           # reached maxit
      phi_j<-phi[j,]
      phi_use_ml_j<-phi_use_ml[j,]
    }
  }
  
  results=list(coefs_j=coefs_j,
               theta_j=theta_j,
               temp_j=temp_j,
               phi_j=phi_j,phi_use_ml_j=phi_use_ml_j)
  return(results)
}




EM<-function(y, k,
             lambda1=0, lambda2=0, tau=0,
             size_factors=rep(1,times=ncol(y)) ,
             norm_y=y,
             true_clusters=NA){
  
  start_time <- Sys.time()
  
  no_cores<-detectCores()-1
  
  n<-ncol(y)
  g<-nrow(y)
  
  # this makes it possible to have y=0 --> adds 0.1 to all y
  y = y+0.1
  
  vect_y<-as.vector(t(y))
  new_y<-rep(vect_y,each=k) # flatten and multiply each count by number of clusters
  gene<-rep(1:g,each=k*n) # gene for each corresponding new_y
  clusts<-matrix(rep(t(diag(k)),times=n*g),byrow=TRUE,ncol=k) # cluster indicators
  
  # EM
  
  # Initial Clustering
  d<-as.dist(1-cor(norm_y, method="spearman"))  ##Spearman correlation distance w/ log transform##
  model<-hclust(d,method="complete")       # hierarchical clustering
  #col<-rep("",times=ncol(y))
  #for(i in 1:length(col)){if(anno$Adeno.Squamous[i]=="adenocarcinoma"){col[i]="red"}else{col[i]="blue"}}
  #heatmap.2(as.matrix(norm_y), Rowv=as.dendrogram(model), Colv=as.dendrogram(model),ColSideColors=col)
  #cls<-cutree(model,k=k)
  #cls<-sample(1:k,n,replace=TRUE) #random initialization
  cls<-cutree(model,k=k)
  
  ########################## SIMULATION ONLY #############################
  if(is.na(true_clusters)==FALSE){
    all_perms=allPerms(1:k)
    all_clusts=list()
    temp_clust<-rep(0,times=n)
    
    for(ii in 1:nrow(all_perms)){
      for(i in 1:n){
        temp_clust[i]<-all_perms[ii,cls[i]]
      }
      all_clusts[[ii]]<-temp_clust
    }
    
    all_clusts[[nrow(all_perms)+1]]<-cls     # contains all permutations of cluster indices
    match_index<-rep(0,times=nrow(all_perms)+1)
    
    for(ii in 1:(nrow(all_perms)+1)){
      match_index[ii]<-mean(true_clusters==all_clusts[[ii]])     # compares each permutation to true --> % of match
    }
    
    cls<-all_clusts[[which.max(match_index)]]
  }
  ########################## SIMULATION ONLY #############################
  
  
  # initialize weights
  wts<-matrix(rep(0,times=k*ncol(y)),nrow=k)
  for(c in 1:k){
    wts[c,]=(cls==c)^2
  }
  
  vect_wts<-rep(as.vector(wts),times=g)
  clust_index<-rep((1:k),times=n*g)
  
  offset=log(size_factors)
  #offset=rep(0,times=n)            # no offsets
  
  dat<-cbind(new_y,clusts,clust_index,gene,vect_wts, rep(rep(offset,each=k),times=g) ) # this is k*g*n rows. cols: count, indicator for cl1, cl2, cl3, genes, wts
  
  colnames(dat)[1]<-c("count")
  colnames(dat)[(k+2):ncol(dat)]<-c("clusts","g","weights","offset")
  
  finalwts<-matrix(rep(0,times=k*ncol(y)),nrow=k)
  
  coefs<-matrix(rep(0,times=g*k),nrow=g)
  pi<-rep(0,times=k)
  theta_list<-list()                # temporary to hold all K x K theta matrices
  
  IRLS_tol = 1E-6                   # Tolerance levels for embedded IRLS and Q fx in EM
  maxit_IRLS = 50
  
  # NB_tol = 1E-5
  # maxit_NB = 100
  
  EM_tol = 1E-6
  maxit_EM = 500
  Q<-rep(0,times=maxit_EM)
  
  lowerK<-0
  
  phi= matrix(0,nrow=g,ncol=k)    # initial gene-specific dispersion parameters for negative binomial
                                  # --> Poisson (phi = 0 = 1/theta)
  #phi=rep(0,times=g)     # gene specific
  
  #glmphi=phi                        # phi's (1/theta) generated from glm
  #init_phi=phi                      # store initialization of phi (1/disp)
  
  temp_list <- list()             # store temp to see progression of IRLS
  
  phi_list <- list()              # store each iteration of phi to see change with each iteration of EM
  
  phi_use_ml = matrix(1,nrow=g,ncol=k)
  
  ########### M / E STEPS #########
  for(a in 1:maxit_EM){
    # M step
    
    dat[,"weights"]<-rep(as.vector(wts),times=g) # update weights column in dat
    
    # IRWLS:
    
    par_X=rep(list(list()),g)
    
    cl<-makeCluster(no_cores)
    i=1                            # cluster complained when this wasn't defined before
    
    clusterExport(cl=cl,varlist=c(ls(),"theta.ml2","soft_thresholding"),envir=environment())
    clusterEvalQ(cl, library("MASS"))
    
    par_X<-parLapply(cl, 1:g, M.step)
    
    stopCluster(cl)
    
    for(j in 1:g){
      coefs[j,] <- par_X[[j]]$coefs_j
      theta_list[[j]] <- par_X[[j]]$theta_j
      temp_list[[j]] <- par_X[[j]]$temp_j
      phi[j,] <- par_X[[j]]$phi_j
      phi_use_ml[j,] <- par_X[[j]]$phi_use_ml_j
    }
    
    phi_list[[a]] <- phi
    
    
    # update on pi_hat
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
    
    
    # E step
    
    # update on weights
    logdenom = apply(log(pi) + l, 2,logsumexpc)
    for(c in 1:k){
      wts[c,]<-exp(log(pi[c])+l[c,]-logdenom)
    }
    
    if(a==maxit_EM){finalwts<-wts}
    # print(pi) # print estimated cluster proportions
    
    
    # if(any(rowSums(wts)==0)){
    #   finalwts<-wts
    #   print(paste("Empty cluster when K =",k,". Choose smaller K"))
    #   lowerK=1
    #   break
    # }
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
  }
  
  
  
  num_warns=length(warnings())
  
  final_clusters<-rep(0,times=n)
  for(i in 1:n){
    final_clusters[i]<-which.max(finalwts[,i])
  }
  
  m<-rep(0,times=g)
  nondiscriminatory=rep(FALSE,times=g)
  
  for(j in 1:g){
    m[j] <- sum(theta_list[[j]][1,]!=0) + 1         # of parameters estimated
    if(m[j]==1){nondiscriminatory[j]=TRUE}
  }
  
  pred.nondiscriminatory<-mean(nondiscriminatory)
  
  
  log_L<-sum(apply(log(pi) + l, 2, logsumexpc))
  
  BIC=-2*log_L+log(n*g)*(sum(m)+(k-1))         # -2log(L) + log(#obs)*(#parameters estimated). minimum = best. g*k: total params, sum(m): total # of discriminatory genes
  if(lowerK==1){BIC=.Machine$integer.max}      # set BIC as max (worst) if K too high
  
  end_time <- Sys.time()
  
  time_elap <- as.numeric(end_time)-as.numeric(start_time)
  
  result<-list(pi=pi,
               coefs=coefs,
               Q=Q[1:a],
               BIC=BIC,
               nondiscriminatory=nondiscriminatory,
               init_clusters=cls,
               final_clusters=final_clusters,
               phi=phi,
               logL=log_L,
               wts=wts,
               time_elap=time_elap)
  return(result)
  
}


