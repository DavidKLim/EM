rm(list=ls())

setwd("/netscr/deelim/out")

source("sim.fx.EM.R")

# run.sim(true_k=2,fold_change=0.3,num_disc=.25,g=100,n=100,method="nb")



run.sim(true_k=c(3,4,5,6,7),method="nb")
