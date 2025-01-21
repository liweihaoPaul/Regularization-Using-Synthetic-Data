# wdbc data
library(readr)
wdbc <- read_csv("wdbc.data", col_names = FALSE)
wdbc<-wdbc[,-1]

Y=as.numeric(wdbc[,1]=="M")
X=as.matrix(wdbc[,-1])
X=X[,1:10]
mode(X)= "numeric"
X=scale(X)

library(glmnet)
library(parallel)
numCores=min(50,detectCores())
library(Rcpp)
library(RcppEigen)
options(rcpp.warnNoExports = FALSE)
sourceCpp("estimate_VE_Eigencpp.cpp")
library(glmtrans)

loocv_best_tau<-function(X,Y,Xstar,Ystar,tau_0_seq){
  p=ncol(X)
  n=length(Y)
  M=length(Ystar)
  ve_error=rep(0,length(tau_0_seq))
  for(tau_0_index in 1:length(tau_0_seq)){
    tau_0=tau_0_seq[tau_0_index]
    fit_cat_tau_0=glmnet(rbind(X,Xstar),c(Y,Ystar),weights =c(rep(1,n),rep(tau_0*n/M,M)) ,
                         family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
    betahat_tau0=as.numeric(fit_cat_tau_0$beta)
    ve_error[tau_0_index]=estimate_VELOOCV_cpp(X,Y,Xstar,Ystar,tau_0,betahat_tau0)
  }
  return(tau_0_seq[which.min(ve_error)])
}


betahat_LOOCV<-function(X_obs,Y_obs,X_anc,Y_anc ){
  tau_0_seq=seq(0,nrow(X_obs),length.out=100)/nrow(X_obs)
  best_tau_0=loocv_best_tau(X_obs,Y_obs,X_anc,Y_anc,tau_0_seq)
  print(best_tau_0)
  fit_cat_besttau_0_loocv=glmnet(rbind(X_obs,X_anc),c(Y_obs,Y_anc),
                                 weights =c(rep(1,length(Y_obs)),rep(best_tau_0*length(Y_obs)/length(Y_anc),length(Y_anc))) ,
                                 family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
  betahat=as.numeric(fit_cat_besttau_0_loocv$beta)
}

# is_separable<-function(X,Y){
#   # if have warning message, return FALSE
#   glm(Y~X-1,family="binomial")$converged
# }

Loocv_synthetic_data<-function(X,Y,NUM_split=20){
  M=400
  p=ncol(X)
  n_target=50
  n_test=nrow(X)-M-n_target
  error_noninformative_loocv=rep(0,NUM_split)
  error_informative_loocv=rep(0,NUM_split)
  error_lasso=rep(0,NUM_split)
  error_ridge=rep(0,NUM_split)
  error_trans=rep(0,NUM_split)
  for(seed in 1:NUM_split){
    print(seed)
    set.seed(123+seed)
    # sample 1:569 without replacement
    index_resch=sample(1:nrow(X))
    
    index_train=index_resch[1:n_target]
    index_test=index_resch[(1+n_target):(n_target+n_test)]
    index_source=index_resch[-(1:(n_target+n_test))]
    X_train=X[index_train,]
    Y_train=Y[index_train]
    X_test=X[index_test,]
    Y_test=Y[index_test]
    X_source=X[index_source,]
    Y_source=Y[index_source]
    
    X_anc_noniform=matrix(rnorm(M*p),M,p)
    Y_anc_noniform=rbinom(M,1,0.5)
    
    betahat_noninformative_loocv=betahat_LOOCV(X_train,Y_train,X_anc_noniform,Y_anc_noniform)
    betahat_informative_loocv=betahat_LOOCV(X_train,Y_train,X_source,Y_source)
    # compute lasso and ridge estimate based on X_train Y_train and cross validation
    betahat_lasso_cv_fit=cv.glmnet(X_train,Y_train,family="binomial",alpha=1,intercept = FALSE,standardize = FALSE)
    betahat_lasso_cv=as.numeric(coef(betahat_lasso_cv_fit,s="lambda.min")[2:(p+1)])
    betahat_ridge_cv_fit=cv.glmnet(X_train,Y_train,family="binomial",alpha=0,intercept = FALSE,standardize = FALSE)
    betahat_ridge_cv=as.numeric(coef(betahat_ridge_cv_fit,s="lambda.min")[2:(p+1)])
    
    Dtrain_traget=list(x=X_train,y=Y_train)
    Dsource=list(list(x=X_source,y=Y_source))
    
    fit.pooled <- glmtrans(target = Dtrain_traget, source = Dsource,
                           family = "binomial", transfer.source.id = "all", intercept = FALSE,standardize = FALSE)
    betahat_transglm=as.numeric(fit.pooled$beta[2:(p+1)])
    
    # investigate the classification error based on X_test Y_test for two betahats
    Yhat_noninformative_loocv=as.numeric(X_test%*%betahat_noninformative_loocv>0)
    Yhat_informative_loocv=as.numeric(X_test%*%betahat_informative_loocv>0)
    Yhat_lasso_cv=as.numeric(X_test%*%betahat_lasso_cv>0)
    Yhat_ridge_cv=as.numeric(X_test%*%betahat_ridge_cv>0)
    Yhat_transglm=as.numeric(X_test%*%betahat_transglm>0)
    error_noninformative_loocv[seed]=sum(Yhat_noninformative_loocv!=Y_test)/length(Y_test)
    error_informative_loocv[seed]=sum(Yhat_informative_loocv!=Y_test)/length(Y_test)
    error_lasso[seed]=sum(Yhat_lasso_cv!=Y_test)/length(Y_test)
    error_ridge[seed]=sum(Yhat_ridge_cv!=Y_test)/length(Y_test)
    error_trans[seed]=sum(Yhat_transglm!=Y_test)/length(Y_test)
  }
  
  
  
  return(list(error_noninformative_loocv,error_informative_loocv ,error_lasso,error_ridge,error_trans ))
}

resultss=Loocv_synthetic_data(X,Y,NUM_split=50 )

mean(resultss[[1]])
mean(resultss[[2]])
mean(resultss[[3]])
mean(resultss[[4]])
mean(resultss[[5]])

sd(resultss[[1]])/sqrt(50)
sd(resultss[[2]])/sqrt(50)
sd(resultss[[3]])/sqrt(50)
sd(resultss[[4]])/sqrt(50)
sd(resultss[[5]])/sqrt(50)

