# check separable
library(MASS)




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
signal_strengthscala = sqrt(n)



for(signal_strengthscalr in c(4.5,6.5,8.5,10.5,12.5,14.5)){
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
    warning_indeicator=NULL
    tryCatch({
      # Attempt to execute the function
      glm(y~X-1,family = "binomial")
      warning_indeicator <<- 0
    }, warning = function(w) {
      # If a warning occurs, set result to 1
      warning_indeicator <<- 1
    })  
      
    return(warning_indeicator)
  }
  print((sapply(1:50,repeat_function)))
}




p = 200
n = 500
s = 40
q = 0.1
M=20*p



delta=n/p
signal_strengthscala = sqrt(n)



for(rho in c(0,0.1,0.2,0.3,0.4)){
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
    beta[signal_index] = sample(c(-1, 1), s, replace = T)*signal_strengthscala
    mu = X %*% beta
    y = rbinom(n, size = 1, prob = f(mu)) 
    warning_indeicator=NULL
    tryCatch({
      # Attempt to execute the function
      glm(y~X-1,family = "binomial")
      warning_indeicator <<- 0
    }, warning = function(w) {
      # If a warning occurs, set result to 1
      warning_indeicator <<- 1
    })  
    
    return(warning_indeicator)
    
  }
  print((sapply(1:50,repeat_function)))
}









