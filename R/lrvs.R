#'@title Perform linear regression variable selection
#'@param X,y data matrix and response
#'@param sigma2 known random error variance
#'@param methods any of: ssvs, bbssl, and debiased_lasso
#'@param alpha significance level, default 0.05
#'@param control control parameters for ssvs, bbssl and debiased_lasso. See ?ssvs, ?BB_SSL, ?SSLasso for control parameters
#'@return a named list of methods. For ssvs, and bbsl: a list of posteriorMean and posterior credible set(a matrix, first column is lower ci, second one is upper ci.) For debiased_lasso: a list of beta_hat ,and confidence_interval.
#'@import BBSSL
#'@export
lrvs = function(X,y,sigma2,methods = c('ssvs','bbssl','debiased_lasso'),
                alpha = 0.05,
                verbose=T,
                ssvs_control = list(),
                bbssl_control = list(),
                debiased_lasso_control = list()){
  n = nrow(X)
  p = ncol(X)
  res = list()
  if('ssvs'%in%methods){
    if(verbose){
      print('Running ssvs')
    }
    ssvs_control = modifyList(ssvs_control_default(),ssvs_control,keep.null = TRUE)
    fit_ssvs = ssvs(X,y,
                    pi0=ssvs_control$pi0,var0=ssvs_control$var0,var1=ssvs_control$var1,
                    sigma2=sigma2,
                    ig_a = ssvs_control$ig_a,
                    ig_b = ssvs_control$ig_b,
                    fix_sigma2 = ssvs_control$fix_sigma2,
                    n_burnin=ssvs_control$n_burnin,
                    n_post=ssvs_control$n_post,
                    printevery=ssvs_control$printevery)
    res$ssvs = list(posteriorMean = c(colMeans(fit_ssvs$beta_draws)),
                    posteriorCS = cbind(apply(fit_ssvs$beta_draws,2,function(z){quantile(z,alpha/2)}),
                                        apply(fit_ssvs$beta_draws,2,function(z){quantile(z,1-alpha/2)})))
  }
  if('bbssl'%in%methods){
    if(verbose){
      print('Running bbssl')
    }
    bbssl_control = modifyList(bbssl_control_default(p),bbssl_control,keep.null = TRUE)
    fit_bb = BB_SSL(y, X, method = bbssl_control$method,
                    lambda=bbssl_control$lambda, NSample=bbssl_control$NSample, a=bbssl_control$a, b=bbssl_control$b,
                    maxiter=bbssl_control$maxiter, eps = bbssl_control$eps, burn.in = bbssl_control$burn.in,
                    length.out = bbssl_control$length.out, discard = bbssl_control$discard,
                    alpha = bbssl_control$alpha, sigma = sqrt(sigma2), initial.beta = bbssl_control$initial.beta,
                    penalty = bbssl_control$penalty, theta=bbssl_control$theta)
    res$bbssl = list(posteriorMean = c(apply(fit_bb$beta,2,mean)),
                     posteriorCS = cbind(apply(fit_bb$beta,2,function(x){quantile(x,alpha/2)}),
                                         apply(fit_bb$beta,2,function(x){quantile(x,1-alpha/2)})))
  }
  if('debiased_lasso'%in%methods){
    if(verbose){
      print('Running debiased lasso')
    }
    debiased_control = modifyList(debiased_control_default(),debiased_lasso_control,keep.null = TRUE)
    fit_debiased = SSLasso(X, y, alpha = alpha, lambda = debiased_control$lambda,
                           mu = debiased_control$mu,
                           intercept = debiased_control$intercept,
                           resol=debiased_control$resol,
                           maxiter=debiased_control$maxiter,
                           threshold = debiased_control$threshold,
                           verbose=debiased_control$verbose)
    res$debiased_lasso = list(beta_hat = fit_debiased$coef,
                              confidence_interval = cbind(fit_debiased$low.lim,fit_debiased$up.lim))
  }
  return(res)
}

debiased_control_default = function(){
  return(list(lambda = NULL,
              mu = NULL,
              intercept = F,
              resol = 1.3,
              maxiter = 50,
              threshold = 1e-2,
              verbose = TRUE))
}

ssvs_control_default = function(){
  return(list(pi0=0.95,
              var0=0.01,
              var1=25,
              fix_sigma2=T,
              ig_a=0.01,
              ig_b=0.01,
              n_burnin=500,
              n_post=2000,
              printevery = 500))
}

bbssl_control_default = function(p){
  return(list(method=3,
              NSample=100,
              maxiter=500,
              lambda = c(7, 0.15),
              a=1,
              b=p,eps = 1e-3, burn.in = FALSE,
              length.out = 50, discard = FALSE, alpha = 3, initial.beta=rep(0,p),
              penalty = "adaptive", theta=0.5))
}
