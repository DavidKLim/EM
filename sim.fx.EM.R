run.sim = function(true_k=c(2:7),fold_change=c(0.1,0.3,0.5,0.7,1),num_nondisc=c(10,25,50,75,90)){
  setwd("/netscr/deelim/out")
  for(i in 1:length(true_k)){for(j in 1:length(fold_change)){for(k in 1:length(num_nondisc)){
      cmd = rep(0, 3)
      cmd[1] = "unlink('.RData') \n source('sim_EM.R') \n"
      cmd[2] = sprintf("X = sim.EM(true.K = %d, fold.change =%f, num.nondisc = %d)\n",
                       true_k[i],fold_change[j],num_nondisc[k])
      out = sprintf("/netscr/deelim/run_sim_%d_%f_%d",true_k[i],fold_change[j],num_nondisc[k])
      out2 = sprintf("%s.out", out)
      cmd[3] = sprintf("save(X, file = '%s')", out2)
      cmdf = paste(cmd, collapse = "")
      write.table(cmdf, file = out, col.names = F, row.names = F, quote = F)
      run = sprintf("bsub -M 4 -q week -o /netscr/deelim/dump R CMD BATCH %s", out)
      system(run)
  }}}
}

collect.sim = function(true_k=c(2:7),fold_change=c(0.1,0.3,0.5,0.7,1),num_nondisc=c(10,25,50,75,90)){
  nsims=length(true_k)*length(fold_change)*length(num_nondisc)
  tab<-matrix(0,nrow=nsims,ncol=10)      # 3 conditions, 7 things to tabulate
  colnames(tab)<-c("true.K","log.fold.change","true.num.disc","K","lambda2","tau","ARI","num.nondisc","sens","false.pos")
  ii=1
  for(i in 1:length(true_k)){for(j in 1:length(fold_change)){for(k in 1:length(num_nondisc)){
      out = sprintf("/netscr/deelim/run_sim_%d_%f_%d",true_k[i],fold_change[j],num_nondisc[k])
      out2 = sprintf("%s.out",out)
      print(out2)
      if(!file.exists(out2)) next
      print(out)
      load(out2)
      tab[ii,1:3]<-c(true_k[i],fold_change[j],num_nondisc[k])
      tab[ii,4]<-X$K
      tab[ii,5]<-X$lambda2
      tab[ii,6]<-X$tau
      tab[ii,7]<-X$ARI
      tab[ii,8]<-X$nondisc
      tab[ii,9]<-X$sens
      tab[ii,10]<-X$falsepos
      ii=ii+1
  }}}
  save(tab,file="final_output.out")
}