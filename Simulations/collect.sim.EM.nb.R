rm(list=ls())

setwd("/netscr/deelim/out")

source("sim.fx.EM.R")

collect.sim(distrib="nb",method="EM",disp="gene")
collect.sim(distrib="nb",method="EM",disp="cluster")