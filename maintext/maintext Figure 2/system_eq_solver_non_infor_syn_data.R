# system of equation solver
# by: Weihao Li
# in this file, we will solve the system of equation in the paper      section: non-informative synthetic data

library(R.utils)
library(sigmoid)
library(MASS)
library(emulator)
library(cubature)
library(mvtnorm) 
library(rootSolve)
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
# non-informative synthetic data ------------------------------------------
# in this section, system of equation contain three unknowns and three equation
# input parameter is given by (delta,m,kappa1,tau_0), output is (sigma_*,gamma_*,alpha_*)
# Two methods are provided, one is using the fixed point iteration, the other is using the Newton Raphson

### method one using Newton Raphson to find the root
Newton_solver_non_infor_syn_data<-function(delta,m,kappa1,tau_0,max_iter=30,tol=1e-03){
  v_old=c(1,1,1)
  v_new=c(0.8,1,1.2)
  cnt=0
  while(norm(v_new-v_old,type="2")>tol & cnt<max_iter){
    v_old=v_new
    JMMres <- withTimeout({
      Jacobian(v_old,delta,m,kappa1,tau_0)
    }, timeout = 300, onTimeout = "warning")
    
    if(is.null(JMMres)){
      return(c(0,0,0))
    }
    GDDres=withTimeout({
      c(NT_f1mc_k20integralkappagiven(v_old,delta,m,kappa1,tau_0),
        NT_f2mc_k20integralkappagiven(v_old,delta,m,kappa1,tau_0),
        NT_f3mc_k20integralkappagiven(v_old,delta,m,kappa1,tau_0))
    }, timeout = 300, onTimeout = "warning")
    if(is.null(GDDres)|is.null(JMMres)){
      return(c(0,0,0))
    }
    v_new=v_old-solve(JMMres,GDDres)
    cnt=cnt+1
    print(v_new)
  }
  if(cnt >= max_iter) {
    warning("Maximum number of iterations reached before convergence!")
  }
  #print(paste0("delta_",delta,"kappa1_",kappa1,"tau_0",tau_0,"finish!"))
  return(v_new)
}

d3rho<-function(x){
  sig_x=sigmoid(x)
  sig_x*(1-sig_x)*(1-2*sig_x)
}


NT_f1mc_k20integralkappagiven<-function(v,delta,m11,kappa1,tau_0,intergal_limit=30){
  sigma=v[1];alpha=v[3];gamma_=v[2]
  hcubature(f1mc_k20integralkappagiven,lowerLimit = c(-intergal_limit,-intergal_limit),
            upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral-sigma^2/2/delta
}
NT_f2mc_k20integralkappagiven<-function(v,delta,m11,kappa1,tau_0,intergal_limit=30){
  sigma=v[1];alpha=v[3];gamma_=v[2]
  hcubature(f2mc_k20integralkappagiven,lowerLimit = c(-intergal_limit,-intergal_limit),
            upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral-1+1/delta-m11
}
NT_f3mc_k20integralkappagiven<-function(v,delta,m11,kappa1,tau_0,intergal_limit=30){
  sigma=v[1];alpha=v[3];gamma_=v[2]
  hcubature(f3mc_k20integralkappagiven,lowerLimit = c(-intergal_limit,-intergal_limit),
            upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral+alpha/2/delta
}

# gradient of f1 f2 and f3 w.r.t to sigma gamma and alpha


GD_NT_f1_integrand_deri_sigma<-function(Z,kappa1,delta,alpha,sigma,gamma_,mm,tau_0){
  gamma_0=tau_0*gamma_/mm
  
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  prox_value_gamma=prox_rho(term1,gamma_)
  prox_value_gamma0=prox_rho(term1,gamma_0)
  
  two_sum=drho(-kappa1*Z[1])*2*(term1-prox_value_gamma)*Z[2]*gamma_*d2rho(prox_value_gamma)/(1+gamma_*d2rho(prox_value_gamma))+
    mm*(term1-prox_value_gamma0)*Z[2]*gamma_0*d2rho(prox_value_gamma0)/(1+gamma_0*d2rho(prox_value_gamma0))-sigma/delta
  
  dmvnorm(Z)*two_sum
}
GD_NT_f1_integrand_deri_gamma<-function(Z,kappa1,delta,alpha,sigma,gamma_,mm,tau_0){
  gamma_0=tau_0*gamma_/mm
  
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  prox_value_gamma=prox_rho(term1,gamma_)
  prox_value_gamma0=prox_rho(term1,gamma_0)
  
  two_sum= -drho(-kappa1*Z[1])  +
    tau_0*(term1-prox_value_gamma0)*Z[2]*drho(prox_value_gamma0)/(1+gamma_0*d2rho(prox_value_gamma0))
  
  dmvnorm(Z)*two_sum
}
GD_NT_f1_integrand_deri_alpha<-function(Z,kappa1,delta,alpha,sigma,gamma_,mm,tau_0){
  gamma_0=tau_0*gamma_/mm
  
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  prox_value_gamma=prox_rho(term1,gamma_)
  prox_value_gamma0=prox_rho(term1,gamma_0)
  
  two_sum=drho(-kappa1*Z[1])*2*(term1-prox_value_gamma)*kappa1*alpha*Z[1]*gamma_*d2rho(prox_value_gamma)/(1+gamma_*d2rho(prox_value_gamma))+
    mm*(term1-prox_value_gamma0)*kappa1*alpha*Z[1]*gamma_0*d2rho(prox_value_gamma0)/(1+gamma_0*d2rho(prox_value_gamma0))-sigma/delta
  
  dmvnorm(Z)*two_sum
}

GD_NT_f2_integrand_deri_sigma<-function(Z,kappa1,delta,alpha,sigma,gamma_,mm,tau_0){
  gamma_0=tau_0*gamma_/mm
  
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  prox_value_gamma=prox_rho(term1,gamma_)
  prox_value_gamma0=prox_rho(term1,gamma_0)
  
  two_sum= -drho(-kappa1*Z[1])*2*gamma_*d3rho(prox_value_gamma)*Z[2]/(1+gamma_*d2rho(prox_value_gamma))^3 -
    mm*gamma_0*d3rho(prox_value_gamma0)*Z[2]/(1+gamma_0*d2rho(prox_value_gamma0))^3
  
  dmvnorm(Z)*two_sum
}

GD_NT_f2_integrand_deri_gamma<-function(Z,kappa1,delta,alpha,sigma,gamma_,mm,tau_0){
  gamma_0=tau_0*gamma_/mm
  
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  prox_value_gamma=prox_rho(term1,gamma_)
  prox_value_gamma0=prox_rho(term1,gamma_0)
  
  two_sum= -drho(-kappa1*Z[1])*(d2rho(prox_value_gamma)-gamma_*d3rho(prox_value_gamma)*drho(prox_value_gamma)/(1+gamma_*d2rho(prox_value_gamma)))/(1+gamma_*d2rho(prox_value_gamma))^2 -
    tau_0*(d2rho(prox_value_gamma0)-gamma_0*d3rho(prox_value_gamma0)*drho(prox_value_gamma0)/(1+gamma_0*d2rho(prox_value_gamma0) ))/(1+gamma_0*d2rho(prox_value_gamma0))^2
  
  dmvnorm(Z)*two_sum
}

GD_NT_f2_integrand_deri_alpha<-function(Z,kappa1,delta,alpha,sigma,gamma_,mm,tau_0){
  gamma_0=tau_0*gamma_/mm
  
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  prox_value_gamma=prox_rho(term1,gamma_)
  prox_value_gamma0=prox_rho(term1,gamma_0)
  
  two_sum= -drho(-kappa1*Z[1])*2*gamma_*d3rho(prox_value_gamma)*alpha*Z[1]/(1+gamma_*d2rho(prox_value_gamma))^3 -
    mm*gamma_0*d3rho(prox_value_gamma0)*alpha*Z[1]/(1+gamma_0*d2rho(prox_value_gamma0))^3
  
  dmvnorm(Z)*two_sum
}

GD_NT_f3_integrand_deri_sigma<-function(Z,kappa1,delta,alpha,sigma,gamma_,mm,tau_0){
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  prox_value_gamma=prox_rho(term1,gamma_)
  
  two_sum= d2rho(-kappa1*Z[1])*Z[2]/(1+gamma_*d2rho(prox_value_gamma))
  
  dmvnorm(Z)*two_sum
}
GD_NT_f3_integrand_deri_gamma<-function(Z,kappa1,delta,alpha,sigma,gamma_,mm,tau_0){
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  prox_value_gamma=prox_rho(term1,gamma_)
  
  two_sum= -d2rho(-kappa1*Z[1])*drho(prox_value_gamma)/(1+gamma_*d2rho(prox_value_gamma))
  
  dmvnorm(Z)*two_sum
}
GD_NT_f3_integrand_deri_alpha<-function(Z,kappa1,delta,alpha,sigma,gamma_,mm,tau_0){
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  prox_value_gamma=prox_rho(term1,gamma_)
  
  two_sum= d2rho(-kappa1*Z[1])*kappa1*alpha*Z[1]/(1+gamma_*d2rho(prox_value_gamma))+1/2/delta
  
  dmvnorm(Z)*two_sum
}

# gradient function is a 3x3 matrix with first row is deriviative of f1 w.r.t sigma gamma alpha
# second row is deriviative of f2 w.r.t sigma gamma alpha
# third row is deriviative of f3 w.r.t sigma gamma alpha
Jacobian<-function(v,delta,m11,kappa1,tau_0,intergal_limit=30){
  sigma=v[1];alpha=v[3];gamma_=v[2]
  Jm=matrix(0,3,3)
  Jm[1,1]=hcubature(GD_NT_f1_integrand_deri_sigma,lowerLimit = c(-intergal_limit,-intergal_limit),
                    upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral
  Jm[1,2]=hcubature(GD_NT_f1_integrand_deri_gamma,lowerLimit = c(-intergal_limit,-intergal_limit),
                    upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral
  Jm[1,3]=hcubature(GD_NT_f1_integrand_deri_alpha,lowerLimit = c(-intergal_limit,-intergal_limit),
                    upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral
  
  Jm[2,1]=hcubature(GD_NT_f2_integrand_deri_sigma,lowerLimit = c(-intergal_limit,-intergal_limit),
                    upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral
  Jm[2,2]=hcubature(GD_NT_f2_integrand_deri_gamma,lowerLimit = c(-intergal_limit,-intergal_limit),
                    upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral
  Jm[2,3]=hcubature(GD_NT_f2_integrand_deri_alpha,lowerLimit = c(-intergal_limit,-intergal_limit),
                    upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral
  
  Jm[3,1]=hcubature(GD_NT_f3_integrand_deri_sigma,lowerLimit = c(-intergal_limit,-intergal_limit),
                    upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral
  Jm[3,2]=hcubature(GD_NT_f3_integrand_deri_gamma,lowerLimit = c(-intergal_limit,-intergal_limit),
                    upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral
  Jm[3,3]=hcubature(GD_NT_f3_integrand_deri_alpha,lowerLimit = c(-intergal_limit,-intergal_limit),
                    upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral
  
  return(Jm)
  
}
# 
# # extra method,  rootSolve for root, very slow
# rootSolve_non_infor_syn_data<-function(delta,m,kappa1,tau_0){
#   v_init=c(1,1,1)
#   have_error22=FALSE
#   tryCatch({root= multiroot(S_minus_Vkappagiven,v_init ,delta=delta,m11=m,kappa1=kappa1,tau_0=tau_0)$root},
#            error=function(e){have_error22<<- TRUE},warning = function(w){have_error22<<- TRUE})
#   if(have_error22==TRUE){
#     v_init=c(kappa1/4,1,1)
#     have_error22=FALSE
#     tryCatch({root= multiroot(S_minus_Vkappagiven,v_init ,delta=delta,m11=m,kappa1=kappa1,tau_0=tau_0)$root},
#              error=function(e){ have_error22<<- TRUE},warning = function(w){have_error22<<- TRUE})
#   }
#   if(have_error22==TRUE){
#     v_init=c(kappa1/5,1,1.2)
#     have_error22=FALSE
#     tryCatch({root= multiroot(S_minus_Vkappagiven,v_init ,delta=delta,m11=m,kappa1=kappa1,tau_0=tau_0)$root},
#              error=function(e){ have_error22<<- TRUE},warning = function(w){have_error22<<- TRUE})
#   }
#   return(root)
# }
# 
# 
f1mc_k20integralkappagiven<-function(Z,kappa1,delta,alpha,sigma,gamma_,mm,tau_0){
  gamma_0=tau_0*gamma_/mm

  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  two_sum=drho(-kappa1*Z[1])*(term1-prox_rho(term1,gamma_))^2+
    mm*drho(0)*(term1-prox_rho(term1,gamma_0))^2

  dmvnorm(Z)*two_sum
}
f2mc_k20integralkappagiven<-function(Z,kappa1,delta,alpha,sigma,gamma_,mm,tau_0){
  gamma_0=tau_0*gamma_/mm

  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  two_sum=2*drho(-kappa1*Z[1])/(1+gamma_*d2rho(prox_rho(term1,gamma_)))+
    2*mm*drho(0)/(1+gamma_0*d2rho(prox_rho(term1,gamma_0)))
  dmvnorm(Z)*two_sum
}
f3mc_k20integralkappagiven<-function(Z,kappa1,delta,alpha,sigma,gamma_,mm,tau_0){
  gamma_0=tau_0*gamma_/mm

  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  two_sum=d2rho(-kappa1*Z[1])*prox_rho(term1,gamma_)
  dmvnorm(Z)*two_sum
}
# 
# S_minus_Vkappagiven<-function(v,delta,m11,kappa1,tau_0,intergal_limit=30){
#   y<- numeric(3)
#   sigma=v[1];alpha=v[3];gamma_=v[2]
#   y[1]=sqrt(2*delta*hcubature(f1mc_k20integralkappagiven,lowerLimit = c(-intergal_limit,-intergal_limit),
#                               upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral)-sigma
#   y[2]=hcubature(f2mc_k20integralkappagiven,lowerLimit = c(-intergal_limit,-intergal_limit),
#                  upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral-1+1/delta-m11
#   y[3]=2*delta*(hcubature(f3mc_k20integralkappagiven,lowerLimit = c(-intergal_limit,-intergal_limit),
#                           upperLimit = c(intergal_limit,intergal_limit),kappa1,delta,alpha,sigma,gamma_,mm=m11,tau_0)$integral)+alpha
#   return(y)
# }

# Second method,  fixed point iteration with monte carlo integration
fixed_point_non_infor_syn_data<-function(delta,m,kappa1,tau_0,maxiter=300,tol=1e0-4,MC_NUM=5000){
  v_old=c(0.8,1,1.2)
  save_20_iterate=matrix(0,20,3)
  for(iter in 1:maxiter){
    Z_matrix=rmvnorm(MC_NUM,rep(0,2))
    v_new = S_vMC(
      v_old,
      delta = delta,
      m = m,
      kappa1 = kappa1,
      tau_0 = tau_0,Z_matrix=Z_matrix,MC_NUM=MC_NUM
    )
    if (norm(v_new - v_old, "2") < tol) {
      break
    }
    v_old = v_new
    save_20_iterate[iter%%21,]=v_new
    #print(v_new)
  }
  # if(iter >= maxiter) {
  #   warning("Maximum number of iterations reached before convergence!")
  # }
  return(apply(save_20_iterate,2,mean))
}
S_vMC<-function(v,delta,m,kappa1,tau_0,Z_matrix,MC_NUM){
  sigma=v[1];gamma_=v[2];alpha=v[3]
  y<- numeric(3)
  y[1]= sqrt(2*delta*mean(sapply(1:MC_NUM,function(ii){
    f1mc(Z_matrix[ii,],kappa1,delta,alpha,sigma,gamma_,m,tau_0)
  })))
  y[2]=gamma_*mean(sapply(1:MC_NUM,function(ii){
    f2mc(Z_matrix[ii,],kappa1,delta,alpha,sigma,gamma_,m,tau_0)
  }))/(1-1/delta+m)
  y[3]= -2*delta*mean(sapply(1:MC_NUM,function(ii){
    f3mc(Z_matrix[ii,],kappa1,delta,alpha,sigma,gamma_,m,tau_0)
  }))
  return(y)
}
f1mc<-function(Z,kappa1,delta,alpha,sigma,gamma_,m,tau_0){
  gamma_0=tau_0*gamma_/m
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  two_sum=drho(-kappa1*Z[1])*(term1-prox_rho(term1,gamma_))^2+
    m*drho(0)*(term1-prox_rho(term1,gamma_0))^2

  two_sum
}

f2mc<-function(Z,kappa1,delta,alpha,sigma,gamma_,m,tau_0){
  gamma_0=tau_0*gamma_/m
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  two_sum=2*drho(-kappa1*Z[1])/(1+gamma_*d2rho(prox_rho(term1,gamma_)))+
    2*m*drho(0)/(1+gamma_0*d2rho(prox_rho(term1,gamma_0)))
  two_sum
}

f3mc<-function(Z,kappa1,delta,alpha,sigma,gamma_,m,tau_0){
  gamma_0=tau_0*gamma_/m
  term1=kappa1*alpha*Z[1]+sigma*Z[2]
  two_sum=d2rho(-kappa1*Z[1])*prox_rho(term1,gamma_)
   
  two_sum
}



# do an experiment to see whether two methods provide same root
# delta=4
# m=20/delta
# kappa1=1
# tau_0=0.5
# Newton_solver_non_infor_syn_data(delta,m,kappa1,tau_0) # take five minutes, 0.7759624 0.9129632 0.6295517
# fixed_point_non_infor_syn_data(delta,m,kappa1,tau_0,MC_NUM=6000) # due to randomness, only closed 0.7718643 0.9098724 0.6317573
# fixed_point_non_infor_syn_data(delta,m,kappa1,tau_0,MC_NUM=5000) # due to randomness, only closed  0.7762872 0.9115252 0.6228384
# 
# tau_0=-0.5 # does not work
# Newton_solver_non_infor_syn_data(delta,m,kappa1,tau_0) 
