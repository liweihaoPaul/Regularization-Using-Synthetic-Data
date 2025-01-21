DSMAP <- function(X, y){
  ## X: design matrix
  ## y: response variable
  n = dim(X)[1]; p = dim(X)[2];M=20*p
  ## Split the data into two halves and run MLE
  sample_index1 <- sample(x = c(1:n), size = 0.5 * n, replace = F)
  sample_index2 <- setdiff(c(1:n), sample_index1)
  # generate synthetic data on two half
  X.syn1=mvrnorm(M, mu = rep(0, p), Sigma = crossprod(X[sample_index1,])/(0.5*n)) 
  X.syn2=mvrnorm(M, mu = rep(0, p), Sigma = crossprod(X[sample_index2,])/(0.5*n)) 
  y.syn1= rbinom(M, 1, rep(0.5,M))
  y.syn2= rbinom(M, 1, rep(0.5,M))
  
  tau_0=p/n*2
  fit1 <- glmnet(rbind(X[sample_index1, ],X.syn1),c(y[sample_index1 ],y.syn1),weights =c(rep(1,n/2),rep(tau_0*n/2/M,M)) ,
                              family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
  fit2 <- glmnet(rbind(X[sample_index2, ],X.syn2),c(y[sample_index2 ],y.syn2),weights =c(rep(1,n/2),rep(tau_0*n/2/M,M)) ,
                 family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
  beta1 <- as.numeric( fit1$beta)
  beta2 <- as.numeric( fit2$beta)
  ## Calculate the variance of each estimator(up to a scaling factor)
  tau1 <- 1/sqrt(diag(solve(t(X[sample_index1, ]) %*% X[sample_index1, ])))
  tau2 <- 1/sqrt(diag(solve(t(X[sample_index2, ]) %*% X[sample_index2, ])))
  M <- beta1*beta2*tau1*tau2
  ## Get the selected feature
  select_index <- analys(M, abs(M), q)
  num_select <- length(select_index)
  return(list(select_index = select_index, num_select = num_select))
}
