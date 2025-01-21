library(glmnet)
library(parallel)
numCores=detectCores()


tau_0=1/4
beta_true_gen="t3"

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

coverage_given_est_kappa1hat<-function(kappa1,delta,tau_0,beta_true_gen="t3",N_rep=100,M_div_p=20,p_seq=c(100 , 400 ,1600)){
  result_matrix_coverage=matrix(0,nc=length(p_seq),nr=N_rep)
  kappa1_seq=seq(0.2,6,0.02)
  load(paste0("result_matrix_kappa1_",kappa1,
              "_delta_",delta,"_tau_0_",
              tau_0,"_beta_true_gen_",beta_true_gen,".RData"))
  load(paste0("TheorySolution_non_infor_syn_data_tau_0_0.25_delta_",delta,".RData"))
  repeat_function_for_eta_square<-function(rep_i){
    cover_result_all_p=rep(0,length(p_seq))
    set.seed(rep_i+100) # for observed data beta
    beta_true_large=generate_beta(p_seq[3],kappa1,beta_true_gen)
    set.seed(rep_i+200) # for observed data X
    X_large=matrix(rnorm(p_seq[3]^2*delta , mean = 0, sd = 1), nc = p_seq[3])
    # set.seed(rep_i+300) # for source data beta and X
    # betas_large=kesee*kappa2/kappa1*beta_true_large+lambda*generate_beta(p_seq[3],kappa1=1,beta_true_gen)
    # X_large=matrix(rnorm(p_seq[3]^2*delta_s , mean = 0, sd = 1), nc = p_seq[3]) # make sure generate different data
    set.seed(rep_i+400) # for synthetic data
    Xstar_large=matrix(rnorm(p_seq[3]^2*M_div_p , mean = 0, sd = 1), nc = p_seq[3])
    for(pindex in 1:length(p_seq)){
      set.seed(rep_i+1000*pindex)
      p=p_seq[pindex]
      n=p*delta
      M=M_div_p*p
      beta_true=beta_true_large[1:p]/sqrt(p)
      
      X = X_large[1:n,1:p]
      Y = rbinom(n, 1, 1 / (1 + exp(-X %*% beta_true)))
      Xstar = Xstar_large[1:M,1:p]
      beta0 = rep(0, p)
      Ystar = rbinom(M, 1, 1 / (1 + exp(-Xstar %*% beta0)))
      fit_cat_tau_0=glmnet(rbind(X,Xstar),c(Y,Ystar),weights =c(rep(1,n),rep(tau_0*n/M,M)) ,
                           family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
      betahat_tau0=as.numeric(fit_cat_tau_0$beta)
      
      beta_true=sqrt(p)*beta_true
      betahat_tau0=sqrt(p)*betahat_tau0
      kappa1_hat=est_kappa1_matrix[rep_i,pindex]
      which_kappa1_index=which.min(abs(kappa1_hat-kappa1_seq))
      TheorySolution_rep_i=TheorySolution[which_kappa1_index,]
      upper_CI=(betahat_tau0+1.96*abs(TheorySolution_rep_i[1]))/abs(TheorySolution_rep_i[3])
      lower_CI=(betahat_tau0-1.96*abs(TheorySolution_rep_i[1]))/abs(TheorySolution_rep_i[3])
      # proportion of CI cover true_beta
      cover_result_all_p[pindex]=mean((beta_true>lower_CI)&(beta_true<upper_CI))
    }
    return(cover_result_all_p)
  }
  
  parallel_eta_list=mclapply(1:N_rep,repeat_function_for_eta_square,mc.cores=numCores)
  for(rep_i in 1:N_rep){
    results_etasquare <- parallel_eta_list[[rep_i]]
    result_matrix_coverage[rep_i,]=results_etasquare
  }
  
  
  # save file with kappa1 delta and tau_0 and beta_true_gen
  save(result_matrix_coverage, file=paste0("result_matrix_coverage_kappa1_",kappa1, "_delta_",delta,"_tau_0_",tau_0,"_beta_true_gen_",beta_true_gen,".RData"))
  
}

for(kappa1 in c(0.5,1,1.5,2)){
  for(delta in c(2,4)){
    coverage_given_est_kappa1hat(kappa1,delta,tau_0,beta_true_gen="t3",N_rep=50,M_div_p=20,p_seq=c(100 , 400 ,1600))
  }
}




# summary

summarycoverage_given_est_kappa1hat<-function(kappa1,delta,tau_0=1/4,beta_true_gen="t3"){
  load(paste0("result_matrix_coverage_kappa1_",kappa1, "_delta_",delta,"_tau_0_",tau_0,"_beta_true_gen_",beta_true_gen,".RData"))
 round( apply(result_matrix_coverage,2,mean),3)[1:2]
  #apply(result_matrix_coverage,2,sd))
  
}
summarycoverage_given_est_kappa1hat(kappa1=0.5,delta=2)
summarycoverage_given_est_kappa1hat(kappa1=1,delta=2)
summarycoverage_given_est_kappa1hat(kappa1=1.5,delta=2)
summarycoverage_given_est_kappa1hat(kappa1=2,delta=2)

summarycoverage_given_est_kappa1hat(kappa1=0.5,delta=4)
summarycoverage_given_est_kappa1hat(kappa1=1,delta=4)
summarycoverage_given_est_kappa1hat(kappa1=1.5,delta=4)
summarycoverage_given_est_kappa1hat(kappa1=2,delta=4)





