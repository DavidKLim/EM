rm(list=ls())

setwd("/netscr/deelim/out")

source("sim.fx.EM.R")

run.sim(distrib="nb",method="EM",disp="gene")
run.sim(distrib="nb",method="EM",disp="cluster")