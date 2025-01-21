# here we will create two dictionary for delta=2 based on tau_0=0.25


rm(list=ls())
source("system_eq_solver_non_infor_syn_data.R")
library(parallel)
numCores <- detectCores()
delta=2
tau_0=1/4
kappa1_seq=seq(0.2,6,0.02)

TheorySolution_each_setting<-function(tau_0,kappa1_seq,delta){
  m=20/delta
  repeat_function<-function(k_index){
    kappa1=kappa1_seq[k_index]
    fixed_point_non_infor_syn_data(delta,m,kappa1,tau_0,MC_NUM=20000)
  }
  TheorySolution<-mclapply(1:length(kappa1_seq),repeat_function,mc.cores=numCores)
  TheorySolution<-do.call(rbind,TheorySolution)
  return(TheorySolution)
}
TheorySolution<-TheorySolution_each_setting(tau_0,kappa1_seq,delta)
save(TheorySolution,file=paste0("TheorySolution_non_infor_syn_data_tau_0_0.25_delta_",delta,".RData"))

correspondence_eta_kappa1=data.frame(kappa1=kappa1_seq,eta_square=0)
for(i in 1:length(kappa1_seq)){
  kappa1=kappa1_seq[i]
  eta_square=TheorySolution[i,1]^2+TheorySolution[i,3]^2*kappa1^2
  correspondence_eta_kappa1[i,2]=eta_square
}
save(correspondence_eta_kappa1,file=
       paste0("correspondence_eta_kappa1_non_infor_syn_data_tau_0_0.25_delta_",delta,".RData"))









