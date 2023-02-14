#'@title stochastic search variable selection
#'@param pi0 prior sparse level
#'@param var0,var1 prior variances for spike, and slab distributions
#'@param fix_sigma2 whether fix random error variance at input
#'@param ig_a,ig_b parameters for inverse gamma prior on sigma2
#'@param n_burnin,n_post burn_in and post draws of MCMC
#'@importFrom mvnfast rmvn
#'@export
ssvs = function(X,y,pi0,var0,var1,sigma2,fix_sigma2=T,ig_a=0.01,ig_b=0.01,n_burnin=1000,n_post=5000,printevery = 10){
  n = nrow(X)
  p = ncol(X)
  XtX = crossprod(X,X)
  Xty = crossprod(X,y)
  beta_draws = matrix(nrow=n_burnin + n_post,ncol=p)
  gamma_draws = matrix(nrow=n_burnin + n_post,ncol=p)
  sigma2_draws = c()
  gamma = rep(0,p)
  if(!fix_sigma2){
    sigma2 = var(y)
  }
  for(i in 1:(n_burnin + n_post)){
    if(i%%printevery==0){
      print(paste('ssvs drawing sample',i))
    }
    d_inv = 1/(gamma*var1 + (1-gamma)*var0)
    A = solve(XtX/sigma2 + diag(d_inv))
    beta = c(rmvn(n=1,mu=A%*%Xty/sigma2,sigma=A))
    if(!fix_sigma2){
      sigma2 = rinvgamma(n/2+ig_a,scale=sum((y-X%*%beta)^2)/2+ig_b)
    }
    for(j in 1:p){
      d1 = (gamma*var1 + (1-gamma)*var0)
      d0 = d1
      d1[j] = var1
      d0[j] = var0
      p1 = log(1-pi0) + sum(dnorm(beta,rep(0,p),sqrt(d1),log = T))
      p0 = log(pi0) + sum(dnorm(beta,rep(0,p),sqrt(d0),log = T))
      p_vec = c(p0,p1)
      p_vec = p_vec-max(p_vec)
      p_vec = exp(p_vec)/sum(exp(p_vec))
      gamma[j] = rbinom(1,1,p_vec[2])
    }
    beta_draws[i,] = beta
    gamma_draws[i,] = gamma
    sigma2_draws[i] = sigma2
  }
  return(list(beta_draws=beta_draws,
              gamma_draws=gamma_draws,
              sigma2_draws=sigma2_draws))
}

rinvgamma = function (n, shape, rate = 1, scale = 1/rate){
  if (missing(rate) && !missing(scale))
    rate <- 1/scale
  1/rgamma(n, shape, rate)
}
