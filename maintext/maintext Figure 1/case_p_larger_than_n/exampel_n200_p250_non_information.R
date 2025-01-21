library(glmnet)
library(parallel)
numCores <- detectCores()
tau_0_seq=seq(0.1,4,0.1)

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

empirical_simulation_non_infor_syn_data<-function(delta,m,kappa1,tau_0,
                                                  design_type="Gaussian",beta_true_gen="normal",p=250,NUM_REPEAT=50){
  n=p*delta
  M=m*n
  tau=tau_0*n
  
  # MSE_store=rep(0,NUM_REPEAT)
  # bias_store=rep(0,NUM_REPEAT)
  # corr_store=rep(0,NUM_REPEAT)
  rep_function<-function(rep_i){
    set.seed(rep_i)
    beta_true=generate_beta(p,kappa1,beta_true_gen)
    X=generate_design_X(n,p,design_type)
    Xstar=matrix(rnorm(p*M,mean=0,sd=sqrt(1/p)),nc=p)
    prop1= 1/(1+exp(-X%*%beta_true))
    y1=rbinom(n,1,prop1)
    prop2= 1/(1+exp(rep(0,M)))
    ystar=rbinom(M,1,prop2)
    wt=c(rep(1,n),rep(tau/M,M))
    fit_lambda_seq<-glmnet(rbind(X,Xstar),c(y1,ystar),weights = wt,
                           family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    beta_hat=as.numeric(fit_lambda_seq$beta )
    # MSE_store[rep_i]=norm(beta_hat-beta_true,type="2")^2/p
    # bias_store[rep_i]=sum(beta_hat*beta_true)/p
    # corr_store[rep_i]=sum(beta_hat*beta_true)/sqrt(sum(beta_hat^2)*sum(beta_true^2))
    return(c(norm(beta_hat-beta_true,type="2")^2/p,sum(beta_hat*beta_true)/sqrt(sum(beta_hat^2)*sum(beta_true^2))))
  }
  parallel_list=mclapply(1:NUM_REPEAT,rep_function,mc.cores=numCores)
  all_rep_matrix=do.call(rbind,parallel_list)
  # return average result and lower and upper bound
  
  return(c(mean(all_rep_matrix[,1]),mean(all_rep_matrix[,1])+1.96/sqrt(NUM_REPEAT)*sd(all_rep_matrix[,1]),mean(all_rep_matrix[,1])-1.96/sqrt(NUM_REPEAT)*sd(all_rep_matrix[,1]),
           mean(all_rep_matrix[,2]),mean(all_rep_matrix[,2])+1.96/sqrt(NUM_REPEAT)*sd(all_rep_matrix[,2]),mean(all_rep_matrix[,2])-1.96/sqrt(NUM_REPEAT)*sd(all_rep_matrix[,2])))
}

for (kappa1 in c(2.5)){
  for (delta in c(0.8)){
    m=20/delta
    design_type="Gaussian"
    beta_true_gen="t3"
    empirical_result_matrix=matrix(0,nrow=length(tau_0_seq),ncol= 6)
    for(tau_0_index in 1:length(tau_0_seq)){
      #print(tau_0_index)
      tau_0=tau_0_seq[tau_0_index]
      empirical_result=empirical_simulation_non_infor_syn_data(delta,m,kappa1,tau_0,
                                                               design_type,beta_true_gen,NUM_REPEAT = 50)
      empirical_result_matrix[tau_0_index,]=empirical_result
    }
    save(empirical_result_matrix,file=paste0("empirical_result_matrix_non_infor_syn_data_delta_",
                                             delta,"_m_",m,"_kappa1_",
                                             kappa1,"_design_type_",
                                             design_type,"_beta_true_gen_",beta_true_gen,".RData"))
  }
}

empirical_result_matrix[,1]
plot(tau_0_seq*200,empirical_result_matrix[,1],type="l",ylim=c(4.5,7.5),xlab=expression(tau),ylab="Square error",col="red")
















