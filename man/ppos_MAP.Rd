% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Imp_fn.R
\name{ppos_MAP}
\alias{ppos_MAP}
\title{Posterior analysis: using MAP prior}
\usage{
ppos_MAP(
  Hy,
  Hy.se,
  tau.dist,
  tau.prior,
  alpha,
  sigma.fs,
  alternative,
  n.burnin = 4000,
  n.iter = 10000
)
}
\arguments{
\item{Hy}{Data set containing treatment effects from historical studies.}

\item{Hy.se}{Data set (historical) containing square of standard error of treatment effects from historical studies.}

\item{tau.dist}{default = NULL}

\item{tau.prior}{Parameters of prior distribution for \code{tau}. default=NULL estimates \code{tau} using DerSimonian and Laird (DL) approach.}

\item{alpha}{Two-sided alpha.}

\item{sigma.fs}{Estimated variance of the treatment effect yet to be observed.}

\item{alternative}{Specify alternative, left tailed = "less" and right tailed = "greater"}

\item{n.burnin}{Number of burn-in iterations in the MCMC.}

\item{n.iter}{Total number of iterations in the MCMC. Post burn-in = n.iter-n.burnin}
}
\value{
summary from JAGS output containing the estimates of all parameters in the Bayesian hierarchical model.
}
\description{
Posterior analysis: using MAP prior
}
