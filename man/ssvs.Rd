% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ssvs.R
\name{ssvs}
\alias{ssvs}
\title{stochastic search variable selection}
\usage{
ssvs(
  X,
  y,
  pi0,
  var0,
  var1,
  sigma2,
  fix_sigma2 = T,
  ig_a = 0.01,
  ig_b = 0.01,
  n_burnin = 1000,
  n_post = 5000,
  printevery = 10
)
}
\arguments{
\item{pi0}{prior sparse level}

\item{var0, var1}{prior variances for spike, and slab distributions}

\item{fix_sigma2}{whether fix random error variance at input}

\item{ig_a, ig_b}{parameters for inverse gamma prior on sigma2}

\item{n_burnin, n_post}{burn_in and post draws of MCMC}
}
\description{
stochastic search variable selection
}
