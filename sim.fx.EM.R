run.sim = function(prefix="",true_k=c(2:7),fold_change=c(0.5,1,2),num_disc=c(.05,.1),g=c(1000,2000),n_per=c(25,50),
                   distrib="nb",method="EM",disp="gene",fixed_parms="F", fixed_coef=6.5,fixed_phi=0.35,
                   ncores=10,nsims=ncores,filt_quant=0.25,filt_method=c("mad","pval","none")){
  setwd("/netscr/deelim/out")

  for(i in 1:length(true_k)){for(j in 1:length(fold_change)){for(k in 1:length(num_disc)){for(l in 1:length(g)){for(m in 1:length(n_per)){
      n = n_per[m]*true_k[i]
      cmd = rep(0, 3)
      cmd[1] = "unlink('.RData') \n source('sim_EM.R') \n"
      if(length(fixed_phi)==1){
        cmd[2] = sprintf("X = sim.EM(true.K = %d, fold.change =%f, num.disc = %f, g = %d,n = %d,
                       distrib='%s', method='%s', filt_quant=%f, filt_method='%s',
                       disp='%s', fixed_parms=%s, fixed_coef=%f, fixed_phi=%f,
                       ncores=%d,nsims=%d,iCluster_compare=T)\n",
                       true_k[i],fold_change[j],num_disc[k],g[l],n,distrib,method,filt_quant,filt_method,disp,fixed_parms,fixed_coef,fixed_phi,ncores,nsims)
      } else{
        cmd[2] = sprintf("X = sim.EM(true.K = %d, fold.change =%f, num.disc = %f, g = %d,n = %d,
                       distrib='%s', method='%s', filt_quant=%f, filt_method='%s',
                       disp='%s', fixed_parms=%s, fixed_coef=%f, fixed_phi=c(%s),
                       ncores=%d,nsims=%d,iCluster_compare=T)\n",
                         true_k[i],fold_change[j],num_disc[k],g[l],n,distrib,method,filt_quant,filt_method,disp,fixed_parms,fixed_coef,paste(fixed_phi,collapse=","),ncores,nsims)
      }
      if(fixed_parms=="F"){
        out = sprintf("/netscr/deelim/run_sim_%s_%s_%s_%s_%d_%f_%f_%d_%d_filt_%s_%f",
                    prefix,distrib,method,disp,true_k[i],fold_change[j],num_disc[k],g[l],n,filt_method,filt_quant)
      } else{
        if(length(fixed_phi)==1){
          out = sprintf("/netscr/deelim/run_sim_%s_%s_%s_%s_%d_%f_%f_%d_%d_fixed_%f_%f_filt_%s_%f",
                      prefix,distrib,method,disp,true_k[i],fold_change[j],num_disc[k],g[l],n,fixed_coef,fixed_phi,filt_method,filt_quant)
        } else{
          out = sprintf("/netscr/deelim/run_sim_%s_%s_%s_%s_%d_%f_%f_%d_%d_fixed_%f_%s_filt_%s_%f",
                        prefix,distrib,method,disp,true_k[i],fold_change[j],num_disc[k],g[l],n,fixed_coef,paste(fixed_phi,collapse="_"),filt_method,filt_quant)
        }
      }
      out2 = sprintf("%s.out", out)
      cmd[3] = sprintf("save(X, file = '%s')", out2)
      cmdf = paste(cmd, collapse = "")
      write.table(cmdf, file = out, col.names = F, row.names = F, quote = F)
      run = sprintf("bsub -M 15 -q week -n %d -R \"span[hosts=1]\" -o /netscr/deelim/dump R CMD BATCH %s", ncores,out)
      system(run)
  }}}}}
}

collect.sim = function(prefix="",true_k=c(2:7),fold_change=c(0.5,1,2),num_disc=c(.05,.1),g=c(1000,2000),n_per=c(25,50),
                       distrib="nb",method="EM",disp="gene",fixed_parms="F", fixed_coef=6.5,fixed_phi=0.35,filt_quant=0.25,filt_method=c("mad","pval","none")){
  nsims=length(true_k)*length(fold_change)*length(num_disc)*length(g)*length(n_per)
  tab<-matrix(0,nrow=nsims,ncol=26)      # 3 conditions, 7 things to tabulate
  colnames(tab)<-c("n","g","log.fold.change","true.K","true.disc","K","disc","lambda2","tau","ARI","sens","false.pos","i_K","i_ARI","i_lambda","pred_acc",
                  "ARI_hc","ARI_med","ARI_EM","ARI_iClust","sil_HC","sil_med","sil_EM","sil_iClust")
  ii=1
  for(i in 1:length(true_k)){for(j in 1:length(fold_change)){for(k in 1:length(num_disc)){for(l in 1:length(g)){for(m in 1:length(n_per)){
      n = n_per[m]*true_k[i]
      if(fixed_parms=="F"){
        out = sprintf("/netscr/deelim/run_sim_%s_%s_%s_%s_%d_%f_%f_%d_%d_filt_%s_%f",
                      prefix,distrib,method,disp,true_k[i],fold_change[j],num_disc[k],g[l],n,filt_method,filt_quant)
      } else{
        if(length(fixed_phi)==1){
          out = sprintf("/netscr/deelim/run_sim_%s_%s_%s_%s_%d_%f_%f_%d_%d_fixed_%f_%f_filt_%s_%f",
                      prefix,distrib,method,disp,true_k[i],fold_change[j],num_disc[k],g[l],n,fixed_coef,fixed_phi,filt_method,filt_quant)
        } else{
          out = sprintf("/netscr/deelim/run_sim_%s_%s_%s_%s_%d_%f_%f_%d_%d_fixed_%f_%s_filt_%s_%f",
                        prefix,distrib,method,disp,true_k[i],fold_change[j],num_disc[k],g[l],n,fixed_coef,paste(fixed_phi,collapse="_"),filt_method,filt_quant)
        }
      }
      out2 = sprintf("%s.out",out)
      print(out2)
      if(!file.exists(out2)) next
      print(out)
      load(out2)
      tab[ii,1:5]<-c(n,g[l],fold_change[j],true_k[i],num_disc[k])
      tab[ii,6]<-X$K
      tab[ii,7]<-X$disc
      tab[ii,8]<-X$lambda2
      tab[ii,9]<-X$tau
      tab[ii,10]<-X$ARI
      tab[ii,11]<-X$sens
      tab[ii,12]<-X$falsepos
      tab[ii,13]<-X$ifinal_k
      tab[ii,14]<-X$imean_ARI
      tab[ii,15]<-X$ifinal_lambda
      tab[ii,16]<-X$mean_pred_acc
      tab[ii,17]<-X$ARI_HC
      tab[ii,18]<-X$ARI_med
      tab[ii,19]<-X$ARI_EM
      tab[ii,20]<-X$ARI_iClust
      tab[ii,21]<-X$sil_HC
      tab[ii,22]<-X$sil_med
      tab[ii,23]<-X$sil_EM
      tab[ii,24]<-X$sil_iClust
      tab[ii,25]<-X$final_K_hc
      tab[ii,26]<-X$final_K_med
      ii=ii+1
  }}}}}
  if(fixed_parms=="F"){
    final_results = sprintf("final_table_%s_%s_%s_filt_%s_%f",distrib,method,disp,filt_method,filt_quant)
  } else{
    if(length(fixed_phi)==1){
      final_results = sprintf("final_table_%s_%s_%s_fixed_%f_%f_filt_%s_%f",distrib,method,disp,fixed_coef,fixed_phi,filt_method,filt_quant)
    } else{
      final_results = sprintf("final_table_%s_%s_%s_fixed_%f_%s_filt_%s_%f",distrib,method,disp,fixed_coef,paste(fixed_phi,collapse="_"),filt_method,filt_quant)
    }
  }
  final_results2 = sprintf("%s.out",final_results)
  save(tab,file=final_results2)
}
