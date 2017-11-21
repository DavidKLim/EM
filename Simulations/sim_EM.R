#setwd("C:/Users/David/Desktop/Research/EM")

init_y<-read.table("init_y.txt")
#init_size_factors<-as.numeric(read.table("init_size_factors.txt")[,1])
#init_norm_y<-read.table("init_norm_y.txt")

init_y<-cbind(init_y,init_y,init_y,init_y,init_y)


library("stats")
library("data.table")
library("DESeq2")
library("mclust")


# # # DESeq analysis
# row_names<-paste("gene",seq(g))
# col_names<-paste("subj",seq(n))
# cts<-as.matrix(y)
# rownames(cts)<-row_names
# colnames(cts)<-col_names
# coldata<-data.frame(matrix(paste("cl",true_clusters,sep=""),nrow=n))
# rownames(coldata)<-colnames(cts)
# colnames(coldata)<-"cluster"
# dds<-DESeqDataSetFromMatrix(countData = cts,
#                             colData = coldata,
#                             design = ~ cluster)
# DESeq_dds<-DESeq(dds)
# size_factors<-sizeFactors(DESeq_dds)
# norm_y<-counts(DESeq_dds,normalized=TRUE)




# Function to simulate data
simulate_data=function(n,k,g,init_pi,b,size_factors,method,phi=matrix(0,nrow=g,ncol=k)){
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
  if(method=="poisson"){
    for(j in 1:g){
      for(i in 1:n){
        y[j,i] = rpois( 1, lambda = size_factors[i]*exp(b[j,cl[i]]))
      }
    }
  } else if(method=="nb"){
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
simulate_data_g=function(n,k,g,init_pi,b,size_factors,method,phi=rep(0,times=g)){
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
  if(method=="poisson"){
    for(j in 1:g){
      for(i in 1:n){
        y[j,i] = rpois( 1, lambda = size_factors[i]*exp(b[j,cl[i]]))
      }
    }
  } else if(method=="nb"){
    for(j in 1:g){
      for(i in 1:n){
        y[j,i] = rnbinom( 1, size = 1/phi[j], mu = size_factors[i]*exp(b[j,cl[i]]))
      }
    }
  }
  result<-list(y=y,z=z)
  return(result)
}


# Function to perform EM on simulated data
sim.EM<-function(true.K, fold.change, num.disc, g, n, method){
  # max n = 100, max #
  
  if(method=="poisson"){
    source("Pan EM.R")
  } else if(method=="nb"){
    source("NB Pan EM par.R")
  } else{
    print("no method input. Defaulting to Poisson")
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
  DESeq_dds<-DESeq(dds)
  init_size_factors<-sizeFactors(DESeq_dds)
  init_norm_y<-counts(DESeq_dds,normalized=TRUE)
  
  # Unpenalized run to find initial cluster estimates based on K=k
  k=true.K
  
  X_init<-EM(y=init_y,k=k,lambda1=0,lambda2=0,tau=0,size_factors=init_size_factors,norm_y=init_norm_y,true_clusters=true_clusters)
  init_coefs<-X_init$coefs              # save init estimates for coefs & pi
  init_pi<-X_init$pi
  init_phi<-X_init$phi
  
  # to prevent error:
   for(j in 1:g){
     for(c in 1:k){
       if(init_coefs[j,c]>12){init_coefs[j,c]=12}
     }
   }
  
  # Mean over clusters, controlled fold change
  
  sim_coefs<-matrix(rep(rowSums(init_coefs)/k,times=k),ncol=k)
  
  fold_change<-fold.change
  nondisc_fold_change<-0         # fixed nondisc fold change
  
  tt<-floor(num.disc*g)
  sim_coefs[1:tt,]<-matrix(rep( fold_change*(c(0:(k-1))+rep((1-k)/2,times=k)) ,times=tt),nrow=tt,byrow=TRUE)+sim_coefs[1:tt,]
  #sim_coefs[(tt+1):g,]<-matrix(rep( nondisc_fold_change*(c(0:(k-1))+rep((1-k)/2,times=k)) ,times=(g-tt)),nrow=(g-tt),byrow=TRUE)+sim_coefs[(tt+1):g,]         # nondisc fold change = 0 so this doesn't get changed

  sim_pi<-rep(1/true.K,times=true.K)
  
  size_factors<-init_size_factors    # use this for all simulations
  
  
  
  
  
  #### SIMULATIONS ####
  
  # Simulations to find K (Order Selection)
  
  sim=100
  choose_k<-rep(0,times=sim)
  
  all_data <- list(list())
  
  for(ii in 1:sim){
    
    # Simulate data based on initial estimates/estimate size factors
    #init_phi=matrix(0.5,nrow=g,ncol=k)
    init_phi=10/exp(sim_coefs)^2
    sim.dat<-simulate_data(n=n,k=true.K,g=g,init_pi=sim_pi,b=sim_coefs,size_factors=size_factors,method=method,phi=init_phi) # same disp param
    #sim.dat<-simulate_data_g(n=n,k=true.K,g=g,init_pi=sim_pi,b=sim_coefs,size_factors=size_factors,method=method,phi=init_phi) # same disp param
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
    
    rowZeroes<-rep(0,times=g)
    for(j in 1:g){
      rowZeroes[j] = mean(y[j,]==0)
    }          # Proportion of zeroes in each gene
    
    idx <- (rowSums(y)>100 & rowZeroes<0.5)      # mark only genes with >100 count total: take out genes with excess 0's or too low count
    
    y <- y[idx,]
    norm_y <- norm_y[idx,]
    true_disc <- true_disc[idx]
    
    all_data[[ii]]<-list(y=y,
                         true_clusters=true_clusters,
                         norm_y=norm_y,
                         true_disc=true_disc)
  }
    
  
  
  
  
  
  
  # Simulation order selection based on unpenalized EM #
  
  K_search=c(2:7)
  list_BIC=matrix(0,nrow=length(K_search),ncol=2)
  list_BIC[,1]=K_search
  for(ii in 1:sim){
    
    print(paste("Iteration",ii,":"))
    y = all_data[[ii]]$y
    true_clusters = all_data[[ii]]$true_clusters
    norm_y = all_data[[ii]]$norm_y
    
    for(aa in 1:nrow(list_BIC)){
      X<-EM(y=y,k=list_BIC[aa,1],lambda1=0,lambda2=0,tau=0,size_factors=size_factors,norm_y=norm_y,true_clusters=true_clusters)  # no penalty
      list_BIC[aa,2]<-X$BIC
      print(list_BIC[aa,])
      print(paste("Cluster agreement:", adjustedRandIndex(X$init_clusters,X$final_clusters)))
      print(paste("Time:",X$time_elap,"seconds"))
    }
    
    choose_k[ii]=list_BIC[which(list_BIC[,2]==min(list_BIC[,2])),1]
  }
  
  tab_k<-table(choose_k)
  
  print(tab_k)
  max_k=as.numeric(rownames(tab_k)[which.max(tab_k)])
  
  
  
  
  
  
  
  
  
  
  # Simulations to do Grid search for tuning parameters
  sim_tune=sim/10
  choose_lambda1<-rep(0,times=sim_tune)
  choose_lambda2<-rep(0,times=sim_tune)
  choose_tau<-rep(0,times=sim_tune)
  
  for(ii in 1:sim_tune){
    
    print(paste("Iteration",ii,"of grid search"))    # track iteration
    
    #simulate data
    y = all_data[[ii]]$y
    true_clusters = all_data[[ii]]$true_clusters
    norm_y = all_data[[ii]]$norm_y
    
    #create matrix for grid search values
    lambda1_search=1
    lambda2_search=seq(from=0.1,to=0.9,by=0.2)
    tau_search=seq(from=0.1,to=0.9,by=0.2)
    
    list_BIC=matrix(0,nrow=length(lambda1_search)*length(lambda2_search)*length(tau_search),ncol=4) # matrix of BIC's: one for each combination of penalty params 
    
    list_BIC[,1]=rep(lambda1_search,each=length(lambda2_search)*length(tau_search))
    list_BIC[,2]=rep(rep(lambda2_search,each=length(tau_search)),times=length(lambda1_search))
    list_BIC[,3]=rep(tau_search,times=length(lambda1_search)*length(lambda2_search))
    
    
    #search for optimal penalty parameters
    for(aa in 1:nrow(list_BIC)){
      X<-EM(y=y,k=k,tau=list_BIC[aa,3],lambda1=list_BIC[aa,1],lambda2=list_BIC[aa,2],size_factors=size_factors,norm_y=norm_y,true_clusters=true_clusters)
      list_BIC[aa,4]<-X$BIC
      print(list_BIC[aa,])
      print(paste("Time:",X$time_elap,"seconds"))
    }
    
    #store optimal penalty parameters
    max_index<-which(list_BIC[,4]==min(list_BIC[,4]))
    choose_tau[ii]<-list_BIC[max_index,3]
    choose_lambda1[ii]<-list_BIC[max_index,1]
    choose_lambda2[ii]<-list_BIC[max_index,2]
    
    print(list_BIC[max_index,])
    
    if(length(max_index)>1){
      warning("more than one max index")
      max_index<-max_index[1]
      choose_tau[ii]<-list_BIC[max_index,3]
      choose_lambda1[ii]<-list_BIC[max_index,1]
      choose_lambda2[ii]<-list_BIC[max_index,2]
      # max_tau<-list_BIC[max_index,3]
      # max_lambda1<-list_BIC[max_index,1]
      # max_lambda2<-list_BIC[max_index,2]
    }
  }
  
  
  #store max penalty parameters based on multiple simulations
  tab_tau<-table(choose_tau)
  tab_lambda1<-table(choose_lambda1)
  tab_lambda2<-table(choose_lambda2)
  print(tab_tau)
  print(tab_lambda1)
  print(tab_lambda2)
  
  max_tau=as.numeric(rownames(tab_tau)[which.max(tab_tau)])[1]
  max_lambda1=as.numeric(rownames(tab_lambda1)[which.max(tab_lambda1)])[1]
  max_lambda2=as.numeric(rownames(tab_lambda2)[which.max(tab_lambda2)])[1]
  
  
  
  
  
  
  
  
  
  
  # Simulations for determining performance based on optimal K and penalty params
  temp_pi<-matrix(rep(0,times=k*sim),nrow=sim)
  #temp_coefs<-list()
  temp_nondisc<-rep(0,times=sim)
  temp_ARI<-rep(0,times=sim)
  temp_sensitivity<-rep(0,times=sim)
  temp_falsepos<-rep(0,times=sim)
  for(ii in 1:sim){
    
    y = all_data[[ii]]$y
    true_clusters = all_data[[ii]]$true_clusters
    norm_y = all_data[[ii]]$norm_y
    true_disc = all_data[[ii]]$true_disc
    
    
    # sometimes errors if wrong K is input
    X<-EM(y=y,k=k,tau=max_tau,lambda1=max_lambda1,lambda2=max_lambda2,size_factors=size_factors,norm_y=norm_y,true_clusters=true_clusters)
    temp_pi[ii,]<-X$pi
    #temp_coefs[[ii]]<-X$coefs
    temp_nondisc[ii]<-mean(X$nondiscriminatory)
    temp_ARI[ii]<-adjustedRandIndex(true_clusters,X$final_clusters)
    if(tt>0){
      temp_sensitivity[ii]<-mean(X$nondiscriminatory[true_disc==TRUE]==FALSE)
    } else {temp_sensitivity[ii]<-NA}
    if(tt<g){
      temp_falsepos[ii]<-mean(X$nondiscriminatory[true_disc==FALSE]==FALSE)
    } else {temp_falsepos[ii]<-NA}
    
    print(paste(temp_nondisc[ii],temp_ARI[ii],temp_sensitivity[ii],temp_falsepos[ii]))
    print(paste("Time:",X$time_elap,"seconds"))
  }
  
  mean_pi<-colSums(temp_pi)/sim
  #mean_coefs<-Reduce('+',temp_coefs)/sim
  mean_nondisc<-mean(temp_nondisc)
  mean_ARI<-mean(temp_ARI)
  mean_sensitivity<-mean(temp_sensitivity)
  mean_falsepos<-mean(temp_falsepos)
  
  
  
  # Store for tabulation:
  results<-list(K=max_k,
                lambda2=max_lambda2,
                tau=max_tau,
                ARI=mean_ARI,
                nondisc=mean_nondisc,
                sens=mean_sensitivity,
                falsepos=mean_falsepos,
                all_data)
  return(results)
}
