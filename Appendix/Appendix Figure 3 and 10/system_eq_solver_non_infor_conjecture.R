# system of equation solver
# by: Weihao Li
# in this file, we will solve the system of equation in the paper      section: conjecture non-informative synthetic data



library(sigmoid)
library(MASS)
library(emulator)
library(cubature)
library(rootSolve)

fixed_point_conjecture_non_infor<-function(delta,kappa1,tau_0,maxiter=300,tol=1e0-4){
  v_old=c(1,1,1)
  for(i in 1:maxiter){
    v_new= S_minus_V_minfty(v_old,delta,kappa1,tau_0)
    if (norm(v_new - v_old, "2") < tol) {
      break
    }
    v_old=v_new
  }
  return(v_new)
}


# rootSolve_conjecture_non_infor<-function(delta,kappa1,tau_0){
#   v_init=c(1,1,1)
#   multiroot(S_minus_V_minfty,v_init ,delta=delta,kappa1=kappa1,tau_0=tau_0)$root
# }

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

f1mc_k20integral_minfty<-function(Z,kappa1,delta,alpha,sigma,gamma_,tau_0){
  
  
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  two_sum=drho(-kappa1*Z[1])*(term1-prox_rho(term1,gamma_))^2
  
  dmvnorm(Z)*two_sum
}
f2mc_k20integral_minfty<-function(Z,kappa1,delta,alpha,sigma,gamma_,tau_0){
  
  
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  two_sum=2*drho(-kappa1*Z[1])/(1+gamma_*d2rho(prox_rho(term1,gamma_)))-
    2*tau_0*gamma_*drho(0)*d2rho(term1)
  dmvnorm(Z)*two_sum
}
f3mc_k20integral_minfty<-function(Z,kappa1,delta,alpha,sigma,gamma_,tau_0){
  
  
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  two_sum=d2rho(-kappa1*Z[1])*prox_rho(term1,gamma_)
  dmvnorm(Z)*two_sum
}
S_minus_V_minfty<-function(v,delta,kappa1,tau_0,intergal_limit=20){
  y<- numeric(3)
  sigma=v[1];alpha=v[3];gamma_=v[2]
  y[1]=sqrt(2*delta*hcubature(f1mc_k20integral_minfty,lowerLimit = c(-intergal_limit,-intergal_limit),
                              upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,tau_0)$integral)
  y[2]=gamma_*hcubature(f2mc_k20integral_minfty,lowerLimit = c(-intergal_limit,-intergal_limit),
                 upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,tau_0)$integral/(1-1/delta)
  y[3]=-2*delta*(hcubature(f3mc_k20integral_minfty,lowerLimit = c(-intergal_limit,-intergal_limit),
                          upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,tau_0)$integral)
  return(y)
}

# test
# fixed_point_conjecture_non_infor(delta=4,kappa1=1,tau_0=0.5) # 0.7291395 0.8972560 0.6194619



