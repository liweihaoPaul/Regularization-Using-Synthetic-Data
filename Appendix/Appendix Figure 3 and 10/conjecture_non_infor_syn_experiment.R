# case: conjecture non informative synthetic data numerical experiments
# do theory first with parallel computing and then empirical version
source("system_eq_solver_non_infor_conjecture.R")
library(parallel)
numCores <- detectCores()
tau_0_seq=seq(0.1,2,0.1)
ConjectureSolution_each_setting<-function(tau_0_seq,kappa1,delta){
  
  repeat_function<-function(tau_0_index){
    
    tau_0=tau_0_seq[tau_0_index]
    roott=fixed_point_conjecture_non_infor(delta,kappa1,tau_0)
    print(paste0("tau_0_index_finish_",tau_0_index))
    return(roott)
  }
  ConjectureSolution<-mclapply(1:length(tau_0_seq),repeat_function,mc.cores=numCores)
  ConjectureSolution<-do.call(rbind,ConjectureSolution)
  return(ConjectureSolution)
}


for (kappa1 in c(0.5,1,1.5,2)){
  for (delta in c(2,4)){
    ConjectureSolution<-ConjectureSolution_each_setting(tau_0_seq,kappa1,delta)
    save(ConjectureSolution,file=paste0("ConjectureSolution_non_infor_syn_data_delta_kappa1_",kappa1,"_delta_",delta,".RData")) 
  }
}

#conduct empirical simulation, write a function to return beta under different setting, Newton method
loglikelihood <- function(X, Y, beta){
  n <- length(Y)
  p <- length(beta)
  loglikelihood <- 0
  for (i in 1:n){
    XiT_beta=sum(X[i,]*beta)
    loglikelihood <- loglikelihood + Y[i] * (XiT_beta) - log(1 + exp(XiT_beta))
  }
  return(loglikelihood)
}
gradient_loglikelihood <- function(X, Y, beta){
  n <- length(Y)
  p <- length(beta)
  gradient <- rep(0, p)
  for (i in 1:n){
    XiT_beta=sum(X[i,]*beta)
    gradient <- gradient + X[i,] * (Y[i] - sigmoid(XiT_beta))
  }
  return(gradient)
}
second_order_derivative_loglikelihood <- function(X, Y, beta){
  n <- length(Y)
  p <- length(beta)
  second_order_derivative <- matrix(0, p, p)
  for (i in 1:n){
    XiT_beta=sum(X[i,]*beta)
    second_order_derivative <- second_order_derivative - sigmoid(XiT_beta) * (1 - sigmoid(XiT_beta)) * (X[i,] %*% t(X[i,]))
  }
  return(second_order_derivative)
}
Expectation_sig_one_min_sig<-function(sigma){
  integral<-function(x){sigmoid(x)*(1-sigmoid(x))*dnorm(x,mean=0,sd=sigma)}
  return(integrate(integral,lower=-Inf,upper=Inf)$value)
}
Expectation_sig_one_min_sig_two_sig_minus_one<-function(sigma){
  integral<-function(x){sigmoid(x)*(1-sigmoid(x))*(2*sigmoid(x)-1)*dnorm(x,mean=0,sd=sigma)}
  return(integrate(integral,lower=-Inf,upper=Inf)$value)
}

# gradient of KL penalty w.r.t beta
gradient_KL_penalty <- function(beta, Var_Xistar){
  sigma=sqrt(norm(beta, type="2")^2*Var_Xistar)
  gradient<-beta*Expectation_sig_one_min_sig(sigma)
  return(-gradient*Var_Xistar)
}
#second order derivative of KL penalty w.r.t beta
second_order_derivative_KL_penalty <- function(beta, Var_Xistar){
  sigma=sqrt(norm(beta, type="2")^2*Var_Xistar)
  p=length(beta)
  second_order_derivative<- diag(p)*Expectation_sig_one_min_sig(sigma)
  +diag(beta)*Expectation_sig_one_min_sig_two_sig_minus_one(sigma)
  return(-second_order_derivative*Var_Xistar)
}
posterior_mode_KL_penalty <- function(X, Y, Var_Xistar,tau){
  n <- length(Y)
  p <- ncol(X)
  beta <- rep(0, p)
  max_iteration <- 30
  tolerance <- 1e-10
  for (i in 1:max_iteration){
    gradient <- gradient_loglikelihood(X, Y, beta) + tau* gradient_KL_penalty(beta, Var_Xistar)
    second_order_derivative <- second_order_derivative_loglikelihood(X, Y, beta) + tau*second_order_derivative_KL_penalty(beta, Var_Xistar)
    beta <- beta - qr.solve(second_order_derivative, gradient)
    if (norm(gradient, type="2") < tolerance){
      break
    }
  }
  if (i == max_iteration){
    print("not converge")
  }
  print(paste0("use iteration ", i, " to find the root"))
  return(beta)
}
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
# write a parallel function to do empirical experiment for each  kappa1, delta
empirical_experiment_each_setting<-function(tau_0,kappa1,delta,p=250,NUM_REPEAT=50,
                                            design_type="Gaussian",beta_true_gen="t3"){
  n=delta*p
  tau=tau_0*n
  repeat_function<-function(rep_i){
    set.seed(rep_i)
    beta_true=generate_beta(p,kappa1,beta_true_gen)
    X=generate_design_X(n,p,design_type)

    prop1= 1/(1+exp(-X%*%beta_true))
    Y=rbinom(n,1,prop1)
    beta_hat=posterior_mode_KL_penalty(X,Y,Var_Xistar=1/p,tau=tau)
    return(c(norm(beta_hat-beta_true,type="2")^2/p,sum(beta_hat*beta_true)/p,
             sum(beta_hat*beta_true)/norm(beta_true,"2")/norm(beta_hat,"2")))
  }
  empirical_solution<-mclapply(1:NUM_REPEAT,repeat_function,mc.cores=numCores)
  empirical_solution<-do.call(rbind,empirical_solution)
  return(apply(empirical_solution,2,mean))
}
for (kappa1 in c(0.5,1,1.5,2)){
#for (kappa1 in c(1,1.5,2)){
  for (delta in c(2,4)){

    design_type="Gaussian"
    beta_true_gen="t3"
    empirical_result_matrix=matrix(0,nrow=length(tau_0_seq),ncol=3)
    for(tau_0_index in 1:length(tau_0_seq)){
      #print(tau_0_index)
      tau_0=tau_0_seq[tau_0_index]
      empirical_result_eachtau=empirical_experiment_each_setting(tau_0,kappa1,delta,p=250,NUM_REPEAT=50,
                                                                 design_type=design_type,beta_true_gen=beta_true_gen)
      empirical_result_matrix[tau_0_index,]=empirical_result_eachtau
    }
    save(empirical_result_matrix,file=paste0("empirical_result_matrix_conjecture_non_infor_syn_data_delta_",
                                             delta,"_kappa1_",
                                             kappa1,"_design_type_",
                                             design_type,"_beta_true_gen_",beta_true_gen,".RData"))
  }
}











