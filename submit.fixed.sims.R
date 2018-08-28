rm(list=ls())

setwd("/netscr/deelim/out")

source("sim.fx.EM.R")

fixed_coefs=c(3.75, 6.5, 7.85) # low to high, based on TCGA BRCA 25, 50, 75 quantiles
fixed_phis=c(0.15, 0.35, 1)

for(i in 1:length(fixed_coefs)){for(j in 1:length(fixed_phis)){
    run.sim(distrib="nb",method="EM",disp="gene",fixed_parms="T",fixed_coef=fixed_coefs[i],fixed_phi=fixed_phis[j],ncores=10,
            filt_method="mad",filt_quant=0.25)
  }}