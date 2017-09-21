run.sim = function(prefix="",true_k=c(2:7),fold_change=c(0.1,0.3,0.5,0.7,1),num_disc=c(10,25,50,75,90)){
  setwd("/netscr/deelim/out")
  for(i in 1:length(true_k)){for(j in 1:length(fold_change)){for(k in 1:length(num_disc)){
      cmd = rep(0, 3)
      cmd[1] = "unlink('.RData') \n source('sim_EM.R') \n"
      cmd[2] = sprintf("X = sim.EM(true.K = %d, fold.change =%f, num.disc = %d)\n",
                       true_k[i],fold_change[j],num_disc[k])
      out = sprintf("/netscr/deelim/run_sim_%s_%d_%f_%d",prefix,true_k[i],fold_change[j],num_disc[k])
      out2 = sprintf("%s.out", out)
      cmd[3] = sprintf("save(X, file = '%s')", out2)
      cmdf = paste(cmd, collapse = "")
      write.table(cmdf, file = out, col.names = F, row.names = F, quote = F)
      run = sprintf("bsub -M 4 -q week -o /netscr/deelim/dump R CMD BATCH %s", out)
      system(run)
  }}}
}

collect.sim = function(prefix="",true_k=c(2:7),fold_change=c(0.1,0.3,0.5,0.7,1),num_disc=c(10,25,50,75,90)){
  nsims=length(true_k)*length(fold_change)*length(num_disc)
  tab<-matrix(0,nrow=nsims,ncol=10)      # 3 conditions, 7 things to tabulate
  colnames(tab)<-c("log.fold.change","true.K","true.disc","K","disc","lambda2","tau","ARI","sens","false.pos")
  ii=1
  for(i in 1:length(true_k)){for(j in 1:length(fold_change)){for(k in 1:length(num_disc)){
      out = sprintf("/netscr/deelim/run_sim_%s_%d_%f_%d",prefix,true_k[i],fold_change[j],num_disc[k])
      out2 = sprintf("%s.out",out)
      print(out2)
      if(!file.exists(out2)) next
      print(out)
      load(out2)
      tab[ii,1:3]<-c(fold_change[j],true_k[i],num_disc[k]/100)
      tab[ii,4]<-X$K
      tab[ii,5]<-(1-X$nondisc)
      tab[ii,6]<-X$lambda2
      tab[ii,7]<-X$tau
      tab[ii,8]<-X$ARI
      tab[ii,9]<-X$sens
      tab[ii,10]<-X$falsepos
      ii=ii+1
  }}}
  save(tab,file="final_output.out")
}
