# library(devtools)
# devtools::install_github("BioinformaticsFMRP/TCGAbiolinks")
# library(TCGAbiolinks)
# query1=GDCquery(project="TCGA-BRCA",
#                 data.category = "Gene expression",
#                 data.type = "Gene expression quantification",
#                 platform = "Illumina HiSeq",  
#                 file.type  = "results",
#                 experimental.strategy = "RNA-Seq",
#                 legacy = TRUE)
# GDCdownload(query1)
# expdat <- GDCprepare(query = query1,
#                      save = TRUE, 
#                      save.filename = "TCGA_BRCA_exp.rda")


library(SummarizedExperiment)
setwd("/Users/limddavid/Documents/Research")

load(file="TCGA_BRCA_exp.rda")

cts <- round(assay(data),round=0)             # default gives raw count
anno <- colData(data)@listData

subtypes <- anno$subtype_PAM50.mRNA
table(subtypes)

idy <- !is.na(subtypes)
known_cts <- cts[,idy]           # 500 subjects

table(subtypes[idy][1:100])

source("comparisons.R")
X2 <- compare(known_cts[1:100,1:100])
