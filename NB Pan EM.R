library(caret)
library(mclust)
library(stats)
library(MASS)
library(fpc)
library(permute)
library(amap)
library(gplots)
library(DESeq2)


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


# theta.ml2() from glm.nb() function from MASS package
# EDIT: use MM if ML diverges
theta.ml2 <-
  function(y, mu, n = sum(weights), weights, limit = 10,
           eps = .Machine$double.eps^0.25,
           trace = FALSE){
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
    while((it <- it + 1) < limit && abs(del) > eps) {
      t0 <- abs(t0)
      del <- score(n, t0, mu, y, weights)/(i <- info(n, t0, mu, y, weights))
      if(del < (-t0) || is.na(del)){
        warning("Theta goes from (+) to (-). Last iteration", it," theta =", signif(t0),". Using method of moments instead")
        t0 = theta.mm(y=y,mu=mu,dfr=n-1,weights=weights)
        break                 # if the delta is changing the sign of t0 from + to -, then break (keep last iteration of t0)
      }
      t0 <- t0 + del
      if(trace) message("theta.ml: iter", it," theta =", signif(t0))
    }
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
    t0
  }



# estimateDispersionsGeneEst() from DESeq2 package
# EDIT: use current phi instead of MM estimates
# source("C:/Users/David/Documents/R/win-library/3.4/DESeq2/unexported.DESeq2.fx.R")
# estimateDispersionsGeneEst2 <- function(object, minDisp=1e-8, kappa_0=1,
#                                        dispTol=1e-6, maxit=100, quiet=TRUE,
#                                        modelMatrix=NULL, niter=1, linearMu=NULL,
#                                        minmu=0.5) {
#   if (!is.null(mcols(object)$dispGeneEst)) {
#     if (!quiet) message("found already estimated gene-wise dispersions, removing these")
#     removeCols <- c("dispGeneEst")
#     mcols(object) <- mcols(object)[,!names(mcols(object)) %in% removeCols,drop=FALSE]
#   }
#   stopifnot(length(minDisp) == 1)
#   stopifnot(length(kappa_0) == 1)
#   stopifnot(length(dispTol) == 1)
#   stopifnot(length(maxit) == 1)
#   if (log(minDisp/10) <= -30) {
#     stop("for computational stability, log(minDisp/10) should be above -30")
#   }
#   
#   # in case the class of the mcols(mcols(object)) are not character
#   object <- sanitizeRowRanges(object)
#   
#   if (is.null(modelMatrix)) {
#     modelMatrix <- getModelMatrix(object) 
#     checkFullRank(modelMatrix)
#     if (nrow(modelMatrix) == ncol(modelMatrix)) {
#       stop("the number of samples and the number of model coefficients are equal,
#            i.e., there are no replicates to estimate the dispersion.
#            use an alternate design formula")
#     }
#     } else {
#       message("using supplied model matrix")
#     }
#   
#   object <- getBaseMeansAndVariances(object)
#   
#   # only continue on the rows with non-zero row mean
#   objectNZ <- object[!mcols(object)$allZero,,drop=FALSE]
#   
#   # this rough dispersion estimate (alpha_hat)
#   # is for estimating mu
#   # and for the initial starting point for line search
#   # first check if model matrix is full rank
#   fullRank <- qr(modelMatrix)$rank == ncol(modelMatrix)
#   alpha_hat <- 1/phi
#   
#   # bound the rough estimated alpha between minDisp and maxDisp for numeric stability
#   maxDisp <- max(10, ncol(object))
#   alpha_hat <- alpha_hat_new <- alpha_init <- pmin(pmax(minDisp, alpha_hat), maxDisp)
#   
#   stopifnot(length(niter) == 1 & niter > 0)
#   
#   # use weights if they are present in assays(object)
#   # (we need this already to decide about linear mu fitting)
#   wlist <- getAndCheckWeights(object, modelMatrix)
#   weights <- wlist$weights
#   useWeights <- wlist$useWeights
#   
#   # use a linear model to estimate the expected counts
#   # if the number of groups according to the model matrix
#   # is equal to the number of columns
#   if (is.null(linearMu)) {
#     modelMatrixGroups <- modelMatrixGroups(modelMatrix)
#     linearMu <- nlevels(modelMatrixGroups) == ncol(modelMatrix)
#     # also check for weights (then can't do linear mu)
#     if (useWeights) {
#       linearMu <- FALSE
#     }
#   }
#   
#   # below, iterate between mean and dispersion estimation (niter) times
#   fitidx <- rep(TRUE,nrow(objectNZ))
#   mu <- matrix(0, nrow=nrow(objectNZ), ncol=ncol(objectNZ))
#   dispIter <- numeric(nrow(objectNZ))
#   # bound the estimated count by 'minmu'
#   # this helps make the fitting more robust,
#   # because 1/mu occurs in the weights for the NB GLM
#   for (iter in seq_len(niter)) {
#     if (!linearMu) {
#       fit <- fitNbinomGLMs(objectNZ[fitidx,,drop=FALSE],
#                            alpha_hat=alpha_hat[fitidx],
#                            modelMatrix=modelMatrix)
#       fitMu <- fit$mu
#     } else {
#       fitMu <- linearModelMuNormalized(objectNZ[fitidx,,drop=FALSE],
#                                        modelMatrix)
#     }
#     fitMu[fitMu < minmu] <- minmu
#     mu[fitidx,] <- fitMu
#     
#     # use of kappa_0 in backtracking search
#     # initial proposal = log(alpha) + kappa_0 * deriv. of log lik. w.r.t. log(alpha)
#     # use log(minDisp/10) to stop if dispersions going to -infinity
#     
#     dispRes <- fitDispWrapper(ySEXP = counts(objectNZ)[fitidx,,drop=FALSE],
#                               xSEXP = modelMatrix,
#                               mu_hatSEXP = fitMu,
#                               log_alphaSEXP = log(alpha_hat)[fitidx],
#                               log_alpha_prior_meanSEXP = log(alpha_hat)[fitidx],
#                               log_alpha_prior_sigmasqSEXP = 1, min_log_alphaSEXP = log(minDisp/10),
#                               kappa_0SEXP = kappa_0, tolSEXP = dispTol,
#                               maxitSEXP = maxit, usePriorSEXP = FALSE,
#                               weightsSEXP = weights, useWeightsSEXP = useWeights)
#     
#     dispIter[fitidx] <- dispRes$iter
#     alpha_hat_new[fitidx] <- pmin(exp(dispRes$log_alpha), maxDisp)
#     # only rerun those rows which moved
#     fitidx <- abs(log(alpha_hat_new) - log(alpha_hat)) > .05
#     alpha_hat <- alpha_hat_new
#     if (sum(fitidx) == 0) break
#   }
#   
#   # dont accept moves if the log posterior did not
#   # increase by more than one millionth,
#   # and set the small estimates to the minimum dispersion
#   dispGeneEst <- alpha_hat
#   if (niter == 1) {
#     noIncrease <- dispRes$last_lp < dispRes$initial_lp + abs(dispRes$initial_lp)/1e6
#     dispGeneEst[which(noIncrease)] <- alpha_init[which(noIncrease)]
#   }
#   dispGeneEstConv <- dispIter < maxit
#   
#   # if lacking convergence from fitDisp() (C++)...
#   refitDisp <- !dispGeneEstConv & dispGeneEst > minDisp*10
#   if (sum(refitDisp) > 0) {
#     
#     dispGrid <- fitDispGridWrapper(y = counts(objectNZ)[refitDisp,,drop=FALSE],
#                                    x = modelMatrix,
#                                    mu = mu[refitDisp,,drop=FALSE],
#                                    logAlphaPriorMean = rep(0,sum(refitDisp)),
#                                    logAlphaPriorSigmaSq = 1, usePrior = FALSE,
#                                    weightsSEXP = weights, useWeightsSEXP = useWeights)
#     dispGeneEst[refitDisp] <- dispGrid
#   }
#   dispGeneEst <- pmin(pmax(dispGeneEst, minDisp), maxDisp)
#   
#   dispDataFrame <- buildDataFrameWithNARows(list(dispGeneEst=dispGeneEst),
#                                             mcols(object)$allZero)
#   mcols(dispDataFrame) <- DataFrame(type=rep("intermediate",ncol(dispDataFrame)),
#                                     description=c("gene-wise estimates of dispersion"))
#   mcols(object) <- cbind(mcols(object), dispDataFrame)
#   assays(object)[["mu"]] <- buildMatrixWithNARows(mu, mcols(object)$allZero)
#   
#   return(object)
# }






EM<-function(y, k,
             lambda1=0, lambda2=0, tau=0,
             size_factors=rep(1,times=ncol(y)) ,
             norm_y=y,
             true_clusters=NA){

n<-ncol(y)
g<-nrow(y)



# Create DESeq Dataset
# cts<-as.matrix(y)
# rownames(cts)<-rownames(y)
# colnames(cts)<-colnames(y)
# coldata<-data.frame(rep(1,times=n))
# rownames(coldata)<-colnames(cts)
# colnames(coldata)<-"dummy"
# dds<-suppressMessages(DESeqDataSetFromMatrix(countData = cts,
#                                              colData = coldata,
#                                              design = ~ 1))
# dds<-suppressMessages(DESeq(dds,quiet=TRUE))




# this makes it possible to have y=0 --> adds 0.1 to all y
y = y+0.1


# Create flattened data matrix
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
  
IRLS_tol = 1E-7                   # Tolerance levels for embedded IRLS and Q fx in EM
maxit_IRLS=100

# NB_tol = 1E-5
# maxit_NB = 100

EM_tol = 1E-6
maxit_EM = 1000
Q<-rep(0,times=maxit_EM)
  
lowerK<-0

phi= matrix(0,nrow=g,ncol=k)    # initial gene-specific dispersion parameters for negative binomial
                                 # --> Poisson (phi = 0 = 1/theta)
#phi=rep(0,times=g)

#glmphi=phi                        # phi's (1/theta) generated from glm
#init_phi=phi                      # store initialization of phi (1/disp)

temp_list <- list()             # store temp to see progression of IRLS
phi_list <- list()              # store each iteration of phi to see change with each iteration of EM

  ########### M / E STEPS #########
    for(a in 1:maxit_EM){
    
  # M step
    
    dat[,"weights"]<-rep(as.vector(wts),times=g) # update weights column in dat
    
      # IRWLS:
    # betaglm<-matrix(0,nrow=g,ncol=k)    # use to compare with glm procedure
    beta<-rep(0,times=k)
    
    for(j in 1:g){
      
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
      #temp<-matrix(0,ncol=(k+1),nrow=maxit_IRLS)      # GENE SPECIFIC
      
      dat_j<-dat[dat[,"g"]==j,]                                  # subset just the j'th gene
      
      for(i in 1:maxit_IRLS){
        eta <- 
          if(a==1 & i==1){
            matrix(rep(beta,times=n),nrow=n,byrow=TRUE)               # first initialization of eta
          }else if(a>1 & i==1){
            matrix(rep(beta,times=n),nrow=n,byrow=TRUE) + offset     # Retrieval of eta for IRLS (prev. beta + offset)
          }else{eta}
        
        
        temp[i,]<-c(beta,phi[j,])              # FOR GENE AND CLUSTER SPECIFIC
        #temp[i,]<-c(beta,phi[j])                # GENE SPECIFIC
        
        for(c in 1:k){
          
          dat_jc<-dat_j[dat_j[,"clusts"]==c,]    # subset j'th gene, c'th cluster

          family=negative.binomial(theta=1/phi[j,c])        # can specify family here (plug in updated phi)
          #family=negative.binomial(theta=1/phi[j])           # GENE SPECIFIC
          
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
          
          #### Maximum Likelihood Estimation (cluster & gene) ####
          phi[j,c] <-
            if(all((dat_jc[dat_jc[,"weights"]==1,"count"]-dat_jc[dat_jc[,"weights"]==1,"count"][1])==0)==FALSE){
              1/theta.ml2(y = dat_jc[,"count"]-0.1,
                        mu = mu[,c],
                        weights = dat_jc[,"weights"],
                        limit=10,trace=FALSE)                # this bypasses error when all counts in cluster are identical or
                                                             # there is just one subject in cluster (this would be 0 disp anyway)
            } else {0}
          ########################################
          
          
          
          #### Method of Moments ####
          # phimm[j,c] <- 1/theta.mm(y = dat_jc[,"count"]-0.1,
          #                        mu=mu[,c],
          #                        dfr = sum(dat_jc[,"weights"])-1,
          #                        weights = dat_jc[,"weights"])
          ###########################
          
          #### Deviance ####
          # phimd[j,c] <- 1/theta.md(y = dat_jc[,"count"]-0.1,
          #                        mu=mu[,c],
          #                        dfr = sum(dat_jc[,"weights"])-1,
          #                        weights = dat_jc[,"weights"])
          ##################
          
          #### CHI SQUARE DAMPENING ####
          # Initialize phi #
          # if(a==1 & i==1){
          #   mu[,c] = exp(eta[,c])
          #   Chi2 = sum(dat_jc[,"weights"]*(dat_jc[,"count"]-mu[,c])^2/mu[,c])
          #   df = sum(dat_jc[,"weights"])-1            # df is number in cluster c, minus 1 (effective df if partial weights)
          #   disp = Chi2/df
          # 
          #   if(disp==0){
          #     warning("Dispersion goes to 0")
          #     disp=1e-10
          #   }
          # 
          #   phi[j,c] = 1/disp
          #   init_phi[j,c] = phi[j,c]
          #   #print(paste("FIRST phi",phi[j,c]))
          # }
          # NB_eta = eta[,c]
          # NB_beta = beta
          # family = negative.binomial(theta=1/phi[j,c])
          # linkinv = family$linkinv
          # mu.eta = family$mu.eta
          # variance = family$variance
          # NB_mu = linkinv(NB_eta)
          # NB_mu.eta.val = mu.eta(NB_eta)
          # good <- (dat_jc[,"weights"]>0) & (NB_mu.eta.val != 0)
          # trans_y <- (NB_eta - offset)[good] + (dat_jc[,"count"][good] - NB_mu[good]) / NB_mu.eta.val[good]
          # w <- sqrt(dat_jc[,"weights"][good]*NB_mu.eta.val[good]^2/variance(NB_mu)[good])
          # NB_beta[c] <-
          #   if(lambda1!=0){
          #     ((lambda1*((sum(NB_beta)-NB_beta[c]) + (sum(theta[c,])-theta[c,c])))  +  ((1/n)*sum(w*trans_y))) / ((lambda1*(k-1)) + (1/n)*sum(w))
          #   } else { sum(w*trans_y)/sum(w) }
          # if(NB_beta[c]<(-100)){
          # warning(paste("Cluster",c,"Gene",j,"goes to -infinity"))
          # NB_beta[c] = -100
          # }
          # if(NB_beta[c]>100){
          #   warning(paste("Cluster",c,"Gene",j,"goes to +infinity"))
          #   NB_beta[c] = 100
          # }
          # 
          # NB_eta<-NB_beta[c] + offset      # add back size factors to eta
          # NB_mu <- linkinv(NB_eta)
          # Chi2 = sum(dat_jc[,"weights"]*(dat_jc[,"count"]-NB_mu)^2/(variance(NB_mu)))
          # disp = Chi2/df            # df doesn't change from before
          # 
          # if(disp==0){
          #   warning("Dispersion goes to 0")
          #   disp=1E-5
          # }
          # phi[j,c] = (disp)*phi[j,c]
          #####################################
          
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
          if(sum((temp[i,]-temp[i-1,])^2)<IRLS_tol){
            coefs[j,]<-beta                                 # reached convergence
            theta_list[[j]]<-theta
            temp_list[[j]]<-temp[1:i,]
            break
          }
        }
        if(i==maxit_IRLS){
          coefs[j,]<-beta
          theta_list[[j]]<-theta  
          temp_list[[j]]<-temp[1:i,]           # reached maxit
        }
        
      }
      
      
      #phi[j] =  1/theta.ml2(y = dat_j[,"count"]-0.1, mu = as.vector(t(mu)), weights = dat_j[,"weights"], limit=10,trace=FALSE)   # GENE SPECIFIC
      
    }
    
    #phi<-1/dispersions(estimateDispersionsGeneEst2(dds))
    
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
        #l[c,i]<-sum(dnbinom(y[,i]-0.1,size=1/phi,mu=exp(coefs[,c] + offset[i]),log=TRUE))         # GENE SPECIFIC
      }    # subtract out 0.1 that was added earlier
    }
    
    
    # store and check Q function
    pt1<-(log(pi)%*%rowSums(wts))
    pt2<-sum(wts*l)
    Q[a]<-pt1+pt2
  
    # break condition for EM
    if(a>10){if(abs(Q[a]-Q[a-10])<EM_tol) {
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
  
result<-list(pi=pi,
             coefs=coefs,
             Q=Q[1:a],
             BIC=BIC,
             nondiscriminatory=nondiscriminatory,
             init_clusters=cls,
             final_clusters=final_clusters,
             phi=phi,
             init_dat=y,
             logL=log_L,
             wts=wts)
return(result)
  
}



