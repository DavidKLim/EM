run.sim = function(prefix="",true_k=c(2,4,6),fold_change=c(1,2),num_disc=c(.05,.1),g=c(2000),n=c(100,200),
                   distrib="nb",method="EM",sim_disp="gene",disp="gene",fixed_parms="F", fixed_coef=8,fixed_phi=0.35,
                   ncores=25,nsims=ncores,filt_quant=0.2,filt_method=c("mad","pval","none"),more_options=""){
  setwd("/pine/scr/d/e/deelim/out")

  for(i in 1:length(true_k)){for(j in 1:length(fold_change)){for(k in 1:length(num_disc)){for(l in 1:length(g)){for(m in 1:length(n)){
      cmd = rep(0, 3)
      cmd[1] = "unlink('.RData') \n source('sim_EM.R') \n"
      if(length(fixed_phi)==1){
        cmd[2] = sprintf("X = sim.EM(true.K = %d, fold.change =%f, num.disc = %f, g = %d,n = %d,
                       distrib='%s', method='%s', filt_quant=%f, filt_method='%s',sim_disp='%s',
                       disp='%s', fixed_parms=%s, fixed_coef=%f, fixed_phi=%f,
                       ncores=%d,nsims=%d,iCluster_compare=T%s)\n",
                       true_k[i],fold_change[j],num_disc[k],g[l],n[m],
                       distrib,method,filt_quant,filt_method,sim_disp,
                       disp,fixed_parms,fixed_coef,fixed_phi,
                       ncores,nsims,more_options)
      } else{
        cmd[2] = sprintf("X = sim.EM(true.K = %d, fold.change =%f, num.disc = %f, g = %d,n = %d,
                       distrib='%s', method='%s', filt_quant=%f, filt_method='%s',sim_disp='%s',
                       disp='%s', fixed_parms=%s, fixed_coef=%f, fixed_phi=c(%s),
                       ncores=%d,nsims=%d,iCluster_compare=T%s)\n",
                         true_k[i],fold_change[j],num_disc[k],g[l],n[m],
                         distrib,method,filt_quant,filt_method,sim_disp,
                         disp,fixed_parms,fixed_coef,paste(fixed_phi,collapse=","),
                         ncores,nsims,more_options)
      }
      if(fixed_parms=="F"){
        fname = sprintf("run_sim_%s_%s_%s_sim%s_%s_%d_%f_%f_%d_%d_filt_%s_%f",
                        prefix,distrib,method,sim_disp,disp,true_k[i],fold_change[j],num_disc[k],g[l],n[m],filt_method,filt_quant)
      } else{
        if(length(fixed_phi)==1){
          fname = sprintf("run_sim_%s_%s_%s_sim%s_%s_%d_%f_%f_%d_%d_fixed_%f_%f_filt_%s_%f",
                          prefix,distrib,method,sim_disp,disp,true_k[i],fold_change[j],num_disc[k],g[l],n[m],fixed_coef,fixed_phi,filt_method,filt_quant)
        } else{
          fname = sprintf("run_sim_%s_%s_%s_sim%s_%s_%d_%f_%f_%d_%d_fixed_%f_%s_filt_%s_%f",
                          prefix,distrib,method,sim_disp,disp,true_k[i],fold_change[j],num_disc[k],g[l],n[m],fixed_coef,paste(fixed_phi,collapse="_"),filt_method,filt_quant)
        }
      }
      out = sprintf("/pine/scr/d/e/deelim/%s",fname)
      out2 = sprintf("%s.out", out)
      cmd[3] = sprintf("save(X, file = '%s')", out2)
      cmdf = paste(cmd, collapse = "")
      write.table(cmdf, file = out, col.names = F, row.names = F, quote = F)
      run = sprintf("sbatch -p general -N 1 -n %d --mem=%dg -t 2- -o /pine/scr/d/e/deelim/dump/%s.log -J %s --wrap='R CMD BATCH %s'", ncores,2*ncores,fname,fname,out)
      Sys.sleep(1)
      system(run)
  }}}}}
}

collect.sim = function(prefix="",true_k=c(2,4,6),fold_change=c(1,2),num_disc=c(.05,.1),g=c(2000),n=c(100,200),
                       distrib="nb",method="EM",sim_disp="gene",disp="gene",fixed_parms="F", fixed_coef=8,fixed_phi=0.35,filt_quant=0.2,filt_method=c("mad","pval","none")){
  nsims=length(true_k)*length(fold_change)*length(num_disc)*length(g)*length(n)
  tab<-matrix(0,nrow=nsims,ncol=38)      # 3 conditions, 7 things to tabulate
  colnames(tab)<-c("n","g","log.fold.change","true.K","true.disc","K","disc","lambda","alpha","ARI",
                   "sens","false.pos","i_K","i_ARI","pred_acc","ARI_hc","ARI_med",
                   "K_HC","K_med","Order_Acc","Filt_sens","Filt_FP",
                   "K_NBMB","ARI_NBMB","K_log_MC","ARI_log_MC","K_vsd_MC","ARI_vsd_MC","K_rld_MC","ARI_rld_MC",
                   "OA_iClust","OA_HC","OA_KM","OA_NBMB","OA_log_MC","OA_vsd_MC","OA_rld_MC","true_order_pred_acc")
  ii=1
  for(i in 1:length(true_k)){for(j in 1:length(fold_change)){for(k in 1:length(num_disc)){for(l in 1:length(g)){for(m in 1:length(n)){
      if(fixed_parms=="F"){
        fname = sprintf("run_sim_%s_%s_%s_sim%s_%s_%d_%f_%f_%d_%d_filt_%s_%f",
                        prefix,distrib,method,sim_disp,disp,true_k[i],fold_change[j],num_disc[k],g[l],n[m],filt_method,filt_quant)
      } else{
        if(length(fixed_phi)==1){
          fname = sprintf("run_sim_%s_%s_%s_sim%s_%s_%d_%f_%f_%d_%d_fixed_%f_%f_filt_%s_%f",
                          prefix,distrib,method,sim_disp,disp,true_k[i],fold_change[j],num_disc[k],g[l],n[m],fixed_coef,fixed_phi,filt_method,filt_quant)
        } else{
          fname = sprintf("run_sim_%s_%s_%s_sim%s_%s_%d_%f_%f_%d_%d_fixed_%f_%s_filt_%s_%f",
                          prefix,distrib,method,sim_disp,disp,true_k[i],fold_change[j],num_disc[k],g[l],n[m],fixed_coef,paste(fixed_phi,collapse="_"),filt_method,filt_quant)
        }
      }
      out = sprintf("/pine/scr/d/e/deelim/%s",fname)
      out2 = sprintf("%s.out",out)
      print(out2)
      if(!file.exists(out2)) next
      print(out)
      load(out2)
      tab[ii,1:5]<-c(n[m],g[l],fold_change[j],true_k[i],num_disc[k])
      tab[ii,6]<-X$K
      tab[ii,7]<-X$disc
      tab[ii,8]<-mean(X$all_lambda)
      tab[ii,9]<-mean(X$all_alpha)
      tab[ii,10]<-X$ARI
      tab[ii,11]<-X$sens
      tab[ii,12]<-X$falsepos
      tab[ii,13]<-X$K_iClust
      tab[ii,14]<-X$ARI_iClust
      tab[ii,15]<-X$mean_pred_acc
      tab[ii,16]<-X$ARI_HC
      tab[ii,17]<-X$ARI_KM
      tab[ii,18]<-X$K_HC
      tab[ii,19]<-X$K_KM
      tab[ii,20]<-mean(X$all_k == true_k[i])
      tab[ii,21]<-X$filt_sens
      tab[ii,22]<-X$filt_falsepos
      tab[ii,23]<-X$K_NBMB
      tab[ii,24]<-X$ARI_NBMB
      tab[ii,25]<-X$K_log_MC
      tab[ii,26]<-X$ARI_log_MC
      tab[ii,27]<-X$K_vsd_MC
      tab[ii,28]<-X$ARI_vsd_MC
      tab[ii,29]<-X$K_rld_MC
      tab[ii,30]<-X$ARI_rld_MC
      Ks_iClust=rep(NA,length(X$all_Ks))
      Ks_HC=rep(NA,length(X$all_Ks))
      Ks_KM=rep(NA,length(X$all_Ks))
      Ks_NBMB=rep(NA,length(X$all_Ks))
      Ks_logMC=rep(NA,length(X$all_Ks))
      Ks_vsdMC=rep(NA,length(X$all_Ks))
      Ks_rldMC=rep(NA,length(X$all_Ks))
      for(id in 1:length(X$all_Ks)){
        Ks_iClust[id] = X$all_Ks[[id]]["iClust"]
        Ks_HC[id] = X$all_Ks[[id]]["HC"]
        Ks_KM[id] = X$all_Ks[[id]]["KM"]
        Ks_NBMB[id] = X$all_Ks[[id]]["NBMB"]
        Ks_logMC[id] = X$all_Ks[[id]]["logMC"]
        Ks_vsdMC[id] = X$all_Ks[[id]]["vsdMC"]
        Ks_rldMC[id] = X$all_Ks[[id]]["rldMC"]
      }
      tab[ii,31]<-mean(Ks_iClust == true_k[i])
      tab[ii,32]<-mean(Ks_HC == true_k[i])
      tab[ii,33]<-mean(Ks_KM == true_k[i])
      tab[ii,34]<-mean(Ks_NBMB == true_k[i])
      tab[ii,35]<-mean(Ks_logMC == true_k[i])
      tab[ii,36]<-mean(Ks_vsdMC == true_k[i])
      tab[ii,37]<-mean(Ks_rldMC == true_k[i])
      if(is.null(X$mean_order_pred_acc)){
        tab[ii,38]=NA
      } else{tab[ii,38]<-X$mean_order_pred_acc}
      ii=ii+1
  }}}}}
  if(fixed_parms=="F"){
    final_results = sprintf("final_table_%s_%s_%s_sim%s_%s_filt_%s_%f",prefix,distrib,method,sim_disp,disp,filt_method,filt_quant)
  } else{
    if(length(fixed_phi)==1){
      final_results = sprintf("final_table_%s_%s_%s_sim%s_%s_fixed_%f_%f_filt_%s_%f",prefix,distrib,method,sim_disp,disp,fixed_coef,fixed_phi,filt_method,filt_quant)
    } else{
      final_results = sprintf("final_table_%s_%s_%s_sim%s_%s_fixed_%f_%s_filt_%s_%f",prefix,distrib,method,sim_disp,disp,fixed_coef,paste(fixed_phi,collapse="_"),filt_method,filt_quant)
    }
  }
  final_results2 = sprintf("%s.out",final_results)
  save(tab,file=final_results2)
}
