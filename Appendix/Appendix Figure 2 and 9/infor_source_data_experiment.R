# case: informative source data numerical experiments
# do theory first with parallel computing and then empirical version
# empirical simulation for informative source function need to be modified
source("system_eq_solver_infor_source_data.R")
library(glmnet)
library(parallel)
numCores <- detectCores()
tau_0_seq=seq(0.1,2,0.1)

TheorySolution_each_setting<-function(tau_0_seq,kappa1,kappa2,delta,kesee){
  m=20/delta
  repeat_function<-function(tau_0_index){
    tau_0=tau_0_seq[tau_0_index]
    fixed_point_infor_source_data(delta,m,kappa1,kappa2,kesee,tau_0)
  }
  TheorySolution<-mclapply(1:length(tau_0_seq),repeat_function,mc.cores=numCores)
  TheorySolution<-do.call(rbind,TheorySolution)
  return(TheorySolution)
}
# run for kappa1=c(0.5,1,1.5,2) ,kappa2=1,kesee=0.9 and delta=c(2,4), save the results with proper name
for (kappa1 in c(0.5,1,1.5,2)){
  for (delta in c(2,4)){
    kappa2 = 1
    kesee=0.9
    TheorySolution<-TheorySolution_each_setting(tau_0_seq,kappa1,kappa2 = 1,delta,kesee)
    save(TheorySolution,file=paste0("TheorySolution_infor_source_data_delta_kappa1_",
                                    kappa1,"_kappa2_",kappa2,"_delta_",delta,"_kesee_",kesee,".RData")) 
  }
}

# conduct empirical simulation, write a function to return beta under different setting
generate_beta<-function(p,kappa1,beta_true_gen="normal"){
  if(beta_true_gen=="normal"){
    return(rnorm(p,sd=kappa1))
  }
  if(beta_true_gen=="uniform"){
    return(runif(p,min=-1,max=1)*kappa1*sqrt(3))
  }
  if(beta_true_gen=="t3"){
    return(rt(p,df=3)*kappa1/sqrt(3))
  }
  if(beta_true_gen=="halfsparse"){
    ind <- sample(c(TRUE, FALSE), p, replace = TRUE)
    v=rnorm(p,sd=kappa1*sqrt(2))
    v[ind]=0
    return(v)
  }
}
generate_design_X<-function(n,p,design_type="Gaussian"){
  if(design_type=="Gaussian"){
    return(matrix(rnorm(p*n,mean=0,sd=sqrt(1/p)),nc=p))
  }
  if(design_type=="t4"){
    return(matrix(rt(p*n,df=4)/sqrt(p*2),nc=p))
  }
  if(design_type=="t3"){
    return(matrix(rt(p*n,df=3)/sqrt(p*3),nc=p))
  }
}
empirical_simulation_infor_source_data<-function(delta,m,kappa1,tau_0,kappa2,kesee,
                                                  design_type="Gaussian",beta_true_gen="normal",p=250,NUM_REPEAT=50){
  n=p*delta
  M=m*n
  tau=tau_0*n

  MSE_store=rep(0,NUM_REPEAT)
  bias_store=rep(0,NUM_REPEAT)
  corr_store=rep(0,NUM_REPEAT)
  for(rep_i in 1:NUM_REPEAT){
    set.seed(rep_i)
    beta_true=generate_beta(p,kappa1,beta_true_gen)
    X=generate_design_X(n,p,design_type)
    Xstar=matrix(rnorm(p*M,mean=0,sd=sqrt(1/p)),nc=p)
    prop1= 1/(1+exp(-X%*%beta_true))
    y1=rbinom(n,1, 1/(1+exp(-X%*%beta_true)))
    lambda=kappa2*(sqrt(1-kesee^2))
    beta0=kesee*kappa2/kappa1*beta_true+lambda*generate_beta(p,kappa1=1,beta_true_gen)
    ystar=rbinom(M,1, 1/(1+exp(-Xstar%*%beta0)))
    wt=c(rep(1,n),rep(tau/M,M))
    fit_lambda_seq<-glmnet(rbind(X,Xstar),c(y1,ystar),weights = wt,
                           family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    beta_hat=as.numeric(fit_lambda_seq$beta )
    MSE_store[rep_i]=norm(beta_hat-beta_true,type="2")^2/p
    bias_store[rep_i]=sum(beta_hat*beta_true)/p
    corr_store[rep_i]=sum(beta_hat*beta_true)/sqrt(sum(beta_hat^2)*sum(beta_true^2))
  }
  # return average result
  return(c(mean(MSE_store),mean(bias_store),mean(corr_store)))
}
# run for kappa1=c(0.5,1,1.5,2) ,kappa2=1,kesee=0.9 and delta=c(2,4), save the results with proper name
for (kappa1 in c(0.5,1,1.5,2)){
  for (delta in c(2,4)){
    m=20/delta
    kappa2 = 1
    kesee=0.9
    design_type="Gaussian"
    beta_true_gen="t3"
    empirical_result_matrix=matrix(0,nrow=length(tau_0_seq),ncol=3)
    for(tau_0_index in 1:length(tau_0_seq)){
      #print(tau_0_index)
      tau_0=tau_0_seq[tau_0_index]
      empirical_result=empirical_simulation_infor_source_data(delta,m,kappa1,tau_0,kappa2,kesee,
                                                              design_type,beta_true_gen,p=250,NUM_REPEAT=50)
      empirical_result_matrix[tau_0_index,]=empirical_result
    }
    save(empirical_result_matrix,file=paste0("empirical_result_matrix_infor_source_data_delta_kappa1_",
                                                 kappa1,"_kappa2_",kappa2,"_delta_",delta,"_kesee_",kesee,".RData"))
  }
}



