ABHqMAP <- function(X, y,X.syn,y.syn){
  n = dim(X)[1]; p = dim(X)[2]
  M=dim(X.syn)[1]
  tau_0=p/n
  fit_cat_tau_0=glmnet(rbind(X,X.syn),c(y,y.syn),weights =c(rep(1,n),rep(tau_0*n/M,M)) ,
                       family = "binomial",alpha=0,lambda = 0,intercept = FALSE,standardize = FALSE)
  betahat_tau0=as.numeric(fit_cat_tau_0$beta)
  var_xbetahat_estimation_result=estimate_eta_square_cpp(X,y,X.syn,y.syn,tau_0,betahat_tau0)
  
  load(paste0("correspondence_eta_kappa1_non_infor_syn_data_tau_0,",tau_0,"_delta_",delta,".RData"))
  load(paste0("TheorySolution_non_infor_syn_data_tau_0,",tau_0,"_delta_",delta,".RData"))
  # load("correspondence_eta_kappa1_non_infor_syn_data_tau_0_0.03_delta_8.33333333333333.RData")
  # load("TheorySolution_non_infor_syn_data_tau_0_0.03_delta_8.33333333333333.RData")
  closed_index=which.min(abs(var_xbetahat_estimation_result-correspondence_eta_kappa1$eta_square ))
  sigma_gamma_alpha=TheorySolution[closed_index,]
  
  # estimate partial correlation
  v_seq=1 / sqrt(diag(solve(t(X) %*% X)))  / sqrt(1 - p/n)
  
  test_statistics=betahat_tau0*sqrt(p/n)/(sigma_gamma_alpha[1]/v_seq)
  pvalues=2*pnorm(-abs(test_statistics))
  # fit <- glm(y ~ X - 1, family = 'binomial', x = TRUE, y = TRUE)
  # adjusted_fit <- adjust_glm(fit, verbose = FALSE, echo = TRUE)
  # pvalues <- summary(adjusted_fit)$coefficients[,4]
  sorted_pvalues = sort(pvalues, decreasing = F, index.return = T)
  if(sum(sorted_pvalues$x <= (1:p)*q/p) > 0){
    BHq_index = max(which(sorted_pvalues$x <= (1:p)*q/p))
    select_index = sorted_pvalues$ix[1:BHq_index]
    num_select <- length(select_index)
    return(list(select_index = select_index, num_select = num_select))
  }else{
    return(list(select_index = NULL, num_select = 0))
  }
}
