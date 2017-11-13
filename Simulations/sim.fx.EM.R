run.sim = function(prefix="",true_k=c(2:7),fold_change=c(0.3,0.5,1),num_disc=c(.05,.15,.25),g=c(100,300),n=c(50,100),method="poisson"){
  setwd("/netscr/deelim/out")
  for(i in 1:length(true_k)){for(j in 1:length(fold_change)){for(k in 1:length(num_disc)){for(l in 1:length(g)){for(m in 1:length(n)){
      cmd = rep(0, 3)
      cmd[1] = "unlink('.RData') \n source('sim_EM.R') \n"
      cmd[2] = sprintf("X = sim.EM(true.K = %d, fold.change =%f, num.disc = %f, g = %d,n = %d, method='%s')\n",
                       true_k[i],fold_change[j],num_disc[k],g[l],n[m],method)
      out = sprintf("/netscr/deelim/run_sim_%s_%d_%f_%f_%d_%d_%s",prefix,true_k[i],fold_change[j],num_disc[k],g[l],n[m],method)
      out2 = sprintf("%s.out", out)
      cmd[3] = sprintf("save(X, file = '%s')", out2)
      cmdf = paste(cmd, collapse = "")
      write.table(cmdf, file = out, col.names = F, row.names = F, quote = F)
      run = sprintf("bsub -M 4 -q week -n 10 -R \"span[hosts=1]\" -o /netscr/deelim/dump R CMD BATCH %s", out)
      system(run)
  }}}}}
}

collect.sim = function(prefix="",true_k=c(2:7),fold_change=c(0.3,0.5,1),num_disc=c(.05,.15,.25),g=c(100,300),n=c(50,100),method="poisson"){
  nsims=length(true_k)*length(fold_change)*length(num_disc)
  tab<-matrix(0,nrow=nsims,ncol=12)      # 3 conditions, 7 things to tabulate
  colnames(tab)<-c("n","g","log.fold.change","true.K","true.disc","K","disc","lambda2","tau","ARI","sens","false.pos")
  ii=1
  for(i in 1:length(true_k)){for(j in 1:length(fold_change)){for(k in 1:length(num_disc)){for(l in 1:length(g)){for(m in 1:length(n)){
      out = sprintf("/netscr/deelim/run_sim_%s_%d_%f_%f_%d_%d_%s",prefix,true_k[i],fold_change[j],num_disc[k],g[l],n[m],method)
      out2 = sprintf("%s.out",out)
      print(out2)
      if(!file.exists(out2)) next
      print(out)
      load(out2)
      tab[ii,1:5]<-c(n[m],g[l],fold_change[j],true_k[i],num_disc[k])
      tab[ii,6]<-X$K
      tab[ii,7]<-(1-X$nondisc)
      tab[ii,8]<-X$lambda2
      tab[ii,9]<-X$tau
      tab[ii,10]<-X$ARI
      tab[ii,11]<-X$sens
      tab[ii,12]<-X$falsepos
      ii=ii+1
  }}}}}
  final_results = sprintf("final_table_%s",method)
  final_results2 = sprintf("%s.out",final_results)
  save(tab,file=final_results2)
}
