source("NB Pan EM par.R")
load("BRCA_env_mad25_pure0.9_3subtypes.RData")

# BRCA_anno: raw annotation file (all samples)
# raw_cts: raw counts (all samples and genes)
# BRCA_all_purity: Purity of all samples
# BRCA_all_subtypes: Annotated subtypes of all samples

# anno: annotation file narrowed down to pre-filtered sample space
# norm_y: Normalized counts (adjusted for sequencing depth) after pre-filtering samples and low-count genes
# y: Raw counts after pre-filtering samples and low-count genes
# BRCA_purity: Purity of pre-filtered samples
# BRCA_subtypes: Annotated subtypes of pre-filtered samples
# mads: MAD scores of all genes
# filt_ids: ID's of genes selected by MAD quantile cutoff (in this case, 25%). Can be changed to a different threshold using 'mads'
# size_factors: Size factors (derived by DESeq2) for all samples (not pre-filtered)
# idy: IDs of pre-filtered samples. size_factors[idy] will give size factors of pre-filtered sample space.

# idx, ids_pure, ids_pure_match, ids_subtype, ids_subtype_match: ID's used to generate data. Can disregard for our purposes


# Change MAD filtering to use just top 5% MAD genes
mad_quant = 5
filt_ids = mads >= quantile(mads,1-mad_quant/100)

# Subset size factors to pre-filtered sample space
size_factors = size_factors[idy]

# Subset y and norm_y to MAD-pre-filtered gene space
y = y[filt_ids,]
norm_y = norm_y[filt_ids,]

# Order selection (OS)
K_max=15       # Max number of clusters to search
K=c(1:K_max)   # search range of order

list_BIC = matrix(NA,nrow=K_max,ncol=2)
list_BIC[,1]=K
colnames(list_BIC) = c("K","BIC")

list_X = list()         # Save outputs for each OS run to use as starting point for penalty parameter search

for(i in 1:nrow(list_BIC)){
  X=EM(y=y, k=K[i],
       lambda=0, alpha=0,
       size_factors=size_factors, norm_y=norm_y,
       purity=rep(1,ncol(y)),             # Adjusting for purity: *(disabled for now)*. Purity all = 1: no adjustments for purity
       offsets=0,                         # Offsets: Can include additional offsets (size_factors already included)
       true_clusters=NA, true_disc=NA,    # for diagnostic tracking, if available
       init_parms=FALSE, init_coefs=matrix(0,nrow=nrow(y),ncol=k), init_phi=matrix(0,nrow=nrow(y),ncol=k), init_cls=NA,   # initial starting point, if available
       n_rinits=5, maxit_inits=15,     # specify # of random initializations to search (recommended: 50), and # of iterations for each search (recommended: 15)
       disp="gene", method="EM",        # disp: "gene" or "cluster" specific dispersion estimates; method: "EM" or "CEM"
       prefix="", dir="NA"              # For creating diagnostic files. Will create directory 'dir' within immediate directory "Diagnostics"
       )
  list_X[[i]]=X
  list_BIC[i,2] = X$BIC
  print(list_BIC[i,])
}

max_k = list_BIC[which.min(list_BIC[,2]),1]
print(paste("OS completed. Optimal order: k =",max_k))

# Penalty Parameters (PP) search
## In clusters, I have each grid point being submitted as a separate job: much faster in massive parallel
## In one script, Have to go through each point line by line: much slower.

lambda_search = seq(0.25,2,0.25)
alpha_search = seq(0,0.5,0.05)      # alpha can't be exactly 1

list_BIC = matrix(NA, nrow=length(lambda_search)*length(alpha_search), ncol=3)
list_BIC[,1] = rep(lambda_search,times=length(alpha_search))
list_BIC[,2] = rep(alpha_search,each=length(lambda_search))

colnames(list_BIC) = c("lambda","alpha","BIC")

for(i in 1:nrow(list_BIC)){
  X = EM(y=y, k=max_k,
         lambda=list_BIC[i,1], alpha=list_BIC[i,2],
         size_factors=size_factors, norm_y=norm_y,
         init_parms=TRUE, init_coefs=list_X[[max_k]]$coefs, init_phi = list_X[[max_k]]$phi, init_cls = list_X[[max_k]]$final_clusters,
         disp="gene", method="EM",
         prefix="", dir="NA")
  
  list_BIC[i,3] = X$BIC
  print(list_BIC[i,])
}

max_ID = which.min(list_BIC[,3])
max_lambda = list_BIC[max_ID,1]
max_alpha = list_BIC[max_ID,2]

print(paste("PP search completed. Optimal parameters: lambda =",max_lambda,", alpha =",max_alpha))


# Final run
# Can choose to save output from PP, but in one script, for the sake of memory, run one more time with optimal parameters:

X = EM(y=y, k=max_k,
       lambda=max_lambda, alpha=max_alpha,
       size_factors=size_factors, norm_y=norm_y,
       init_parms=TRUE, init_coefs=list_X[[max_k]]$coefs, init_phi = list_X[[max_k]]$phi, init_cls = list_X[[max_k]]$final_clusters,
       disp="gene", method="EM",
       prefix="", dir="NA")

library(mclust)
adjustedRandIndex(X$final_clusters,BRCA_subtypes)