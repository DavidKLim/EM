rm(list=ls())

setwd("/netscr/deelim/out")

source("sim.fx.EM.R")

run.sim(method="nb",true_k=2,fold_change=0.3,num_disc=75)