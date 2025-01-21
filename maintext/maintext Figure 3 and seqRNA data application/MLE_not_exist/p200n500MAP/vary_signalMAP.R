rm(list = ls())
library(MASS)
library(glmnet)
library(parallel)
num_cores <- detectCores()
library(Rcpp)
library(RcppEigen)
sourceCpp("estimate_eta_square_Eigencpp.cpp")

### source code
source("ABHqMAP.R")
source("ABYqMAP.R")
source("DSMAP.R")
source("MDSMAP.R")
source("utilsMAP.R")


f <- function(x){
  exp(x)/(1 + exp(x))
}

f <- function(x){
  exp(x)/(1 + exp(x))
}



p = 200
n = 500
s = 40
q = 0.1
M=20*p


delta=n/p
rho=0.2

sig_sequence=sqrt(n)*c(1,1.25,1.5,1.75,2)

for(signal_strengthscalr in sig_sequence){
  repeat_function<-function(rep_i){
    set.seed(rep_i)
    
    Sigma = matrix(0, nrow = p, ncol = p)
    for(i in 1:p){
      for(j in 1:p){
        Sigma[i,j] = rho^abs(i-j)
      }
    }
    X = mvrnorm(n, mu = rep(0, p), Sigma = Sigma)*1/sqrt(n)
    signal_index = sample(1:p, size = s, replace = F)
    beta = numeric(p)
    beta[signal_index] = sample(c(-1, 1), s, replace = T)*signal_strengthscalr
    mu = X %*% beta
    y = rbinom(n, size = 1, prob = f(mu)) 
    
    # generate synthetic data
    sigma_hat=crossprod(X)/n
    X.syn = mvrnorm(M, mu = rep(0, p), Sigma = sigma_hat)
    y.syn = rbinom(M, 1, rep(0.5,M))
    
    
    MDS_time = system.time(MDS_result <- MDSMAP(X, y, num_split = 50))[3]
    MDS_result <- fdp_power(MDS_result$select_index, signal_index)
    MDS_fdp <- MDS_result$fdp
    MDS_power <- MDS_result$power
    
    ABHq_time = system.time(ABHq_result <- ABHqMAP(X, y,X.syn,y.syn))[3]
    ABHq_result <- fdp_power(ABHq_result$select_index, signal_index)
    ABHq_fdp <- ABHq_result$fdp
    ABHq_power <- ABHq_result$power
    
    ABYq_time = system.time(ABYq_result <- ABYqMAP(X, y,X.syn,y.syn))[3]
    ABYq_result <- fdp_power(ABYq_result$select_index, signal_index)
    ABYq_fdp <- ABYq_result$fdp
    ABYq_power <- ABYq_result$power
    
    
    ### save data
    data_save <- list(#DS_fdp   = DS_fdp,   DS_power = DS_power, DS_time = DS_time,
      repindex = rep_i, rho = rho, p = p, n = n, s = s, q = q,
      MDS_fdp  = MDS_fdp,  MDS_power = MDS_power, MDS_time = MDS_time,
      ABHq_fdp = ABHq_fdp, ABHq_power = ABHq_power, ABHq_time = ABHq_time,
      ABYq_fdp = ABYq_fdp, ABYq_power = ABYq_power, ABYq_time = ABYq_time)
    return(data_save)
  }
  data_save_list=mclapply(1:50, repeat_function, mc.cores = num_cores)
  save(data_save_list, file = paste0('May3_p',p,'MAPvary_signal_',signal_strengthscalr,'s',s,'cor',rho,'.RData'))
}



# summary
for(signal_strengthscalr in sig_sequence){
  load(paste0('May3_p',p,'MAPvary_signal_',signal_strengthscalr,'s',s,'cor',rho,'.RData'))
  result_matrix=matrix(unlist(data_save_list),nc=15,byrow = T)
  fdr_col_index=c(7,10,13 )
  power_col_index=c(8,11,14 )
  print(rho)

  print( apply(result_matrix[,fdr_col_index],2,mean))
  print( apply(result_matrix[,power_col_index],2,mean))
}
















