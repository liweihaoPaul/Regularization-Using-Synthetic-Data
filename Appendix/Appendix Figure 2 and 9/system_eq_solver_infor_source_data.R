# system of equation solver
# by: Weihao Li
# in this file, we will solve the system of equation in the paper      section: informative source data
# there are four unknowns and four equations, we will use fixed point iteration to solve the system of equation monte carlo integration

library(mvtnorm) 
library(R.utils)
library(sigmoid)
library(MASS)
library(emulator)
rho<-function(x){
  if(x>10) return(10.00005)
  log(1+exp(x))
}
drho<-function(x){
  sigmoid(x)
}
d2rho<-function(x){
  sigmoid(x)*(1-sigmoid(x))
}
prox_rho<-function(v,t){          # argmin_x t*pho(x)+|v-x|^2/2
  
  return(optim(0,function(x){(x-v)^2/2+t*rho(x)},method = "Brent",lower = -1000000000,upper = 1000000000)$par)
  
  newf<-function(xx){
    xx+t*drho(xx)-v
  }
  uniroot(newf,c(-10000,10000))$root
}

fixed_point_infor_source_data<-function(delta,m,kappa1,kappa2,kesee,tau_0,maxiter=200,tol=1e0-4,MC_NUM=20000){
  v_old=c(0.8,0.8,0.8,0.8)
  save_20_iterate=matrix(0,20,4)
  for(iter in 1:maxiter){
    
    v_new = S_vMC_infor(
      v_old,
      delta = delta,
      m = m,
      kappa1 = kappa1,kappa2 = kappa2,kesee = kesee,
      tau_0 = tau_0,MC_NUM=MC_NUM
    )
    if (norm(v_new - v_old, "2") < tol) {
      break
    }
    v_old = v_new
    #print(v_new)
  }
  MC_NUM=50000
  for(iter2 in 1:20){
    v_old = v_new
    save_20_iterate[iter2,] = S_vMC_infor(
      v_old,
      delta = delta,
      m = m,
      kappa1 = kappa1,kappa2 = kappa2,kesee = kesee,
      tau_0 = tau_0,MC_NUM=MC_NUM
    )
    #print(save_20_iterate[iter2,])
  }
  # if(iter >= maxiter) {
  #   warning("Maximum number of iterations reached before convergence!")
  # }
  return(apply(save_20_iterate,2,mean))
}
S_vMC_infor<-function(v,delta,m,kappa1,kappa2,kesee,tau_0,MC_NUM){
  sigma=v[1];gamma_=v[2];alpha1=v[3];alpha2=v[4]
  y<- numeric(4)
  Z_matrix=rmvnorm(MC_NUM,rep(0,3))
  y[1]= sqrt(2*delta*mean(sapply(1:MC_NUM,function(ii){
    f1MC(Z_matrix[ii,],kappa1,kappa2,delta,alpha1,alpha2,sigma,gamma_,m,tau_0,kesee)
  })))
  y[2]=gamma_*mean(sapply(1:MC_NUM,function(ii){
    f2MC(Z_matrix[ii,],kappa1,kappa2,delta,alpha1,alpha2,sigma,gamma_,m,tau_0,kesee)
  }))/(1-1/delta+m)
  y[3]= -2*delta*mean(sapply(1:MC_NUM,function(ii){
    f3MC(Z_matrix[ii,],kappa1,kappa2,delta,alpha1,alpha2,sigma,gamma_,m,tau_0,kesee)
  }))
  y[4]=-2*delta*(mean(sapply(1:MC_NUM,function(ii){
    f4MC(Z_matrix[ii,],kappa1,kappa2,delta,alpha1,alpha2,sigma,gamma_,m,tau_0,kesee)
  })))
  
  return(y)
}
f1MC<-function(Z,kappa1,kappa2,delta,alpha1,alpha2,sigma,gamma_,m,tau_0,kesee){
  
  gamma_0=tau_0*gamma_/m
  
  term1=kappa1*alpha1*Z[1]+kappa2*alpha2*Z[2]+sigma*Z[3]
  two_sum=drho(-kappa1*Z[1])*(term1-prox_rho(term1,gamma_))^2+
    m*drho(-kappa2*kesee*Z[1]-kappa2*sqrt(1-kesee^2)*Z[2])*(term1-prox_rho(term1,gamma_0))^2
  
  two_sum
}

f2MC<-function(Z,kappa1,kappa2,delta,alpha1,alpha2,sigma,gamma_,m,tau_0,kesee){
  gamma_0=tau_0*gamma_/m
  term1=kappa1*alpha1*Z[1]+kappa2*alpha2*Z[2]+sigma*Z[3]
  two_sum=2*drho(-kappa1*Z[1])/(1+gamma_*d2rho(prox_rho(term1,gamma_)))+
    2*m*drho(-kappa2*kesee*Z[1]-kappa2*sqrt(1-kesee^2)*Z[2])/(1+gamma_0*d2rho(prox_rho(term1,gamma_0)))
  two_sum
}

f3MC<-function(Z,kappa1,kappa2,delta,alpha1,alpha2,sigma,gamma_,m,tau_0,kesee){
  gamma_0=tau_0*gamma_/m
  term1=kappa1*alpha1*Z[1]+kappa2*alpha2*Z[2]+sigma*Z[3]
  two_sum=d2rho(-kappa1*Z[1])*prox_rho(term1,gamma_)+
    m*kesee*kappa2/kappa1*d2rho(-kappa2*kesee*Z[1]-kappa2*sqrt(1-kesee^2)*Z[2])*prox_rho(term1,gamma_0)
  two_sum
}

f4MC<-function(Z,kappa1,kappa2,delta,alpha1,alpha2,sigma,gamma_,m,tau_0,kesee){
  gamma_0=tau_0*gamma_/m
  term1=kappa1*alpha1*Z[1]+kappa2*alpha2*Z[2]+sigma*Z[3]
  one_sum=m*sqrt(1-kesee^2)*d2rho(-kappa2*kesee*Z[1]-kappa2*sqrt(1-kesee^2)*Z[2])*prox_rho(term1,gamma_0)
  one_sum
}

# test 
#fixed_point_infor_source_data(delta=4,m=5,kappa1=1,kappa2=1,kesee=0.9,tau_0=0.5)



