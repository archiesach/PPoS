

#' Posterior analysis: Fits bivariate Bayesian hierarchical model to obtain the samples from posterior predictive distribution
#' of \eqn{T_f} to obtain the predictive probability of success (PPoS) given the critical region R from \code{decision_boundary}
#'
#' @param Hy Data set (historical) containing both final and surrogate endpoint from historical studies. A matrix of dimension K x 2 with first and second column containing final and surrogate endpoint treatment effect, respectively.
#' @param Hy.se Data set (historical) containing square of standard error of both final and surrogate endpoint from historical studies. A matrix of dimension K x 2 with first and second column containing final and surrogate endpoint variance (study specific), respectively.
#' @param rho.hat Estimate of study specific correlation coefficient.
#' @param tau.dist Type of prior distribution for between-study variance of both final and surrogate endpoint; supported priors are \code{HalfNormal} (default), \code{Uniform}, \code{InvGamma}, \code{LogNormal}, and \code{HalfCauchy}.
#' @param tau.prior Parameters of prior distribution for \code{tau}
#' @param alpha Two-sided alpha.
#' @param X.cs Observed treatment effect of surrogate endpoint in ongoing study.
#' @param sigma.fs Estimated variance of the treatment effect yet to be observed.
#' @param sigma.cs Variance of the treatment effect of surrogate endpoint observed in the ongoing trial.
#' @param alternative Specify alternative, left tailed = less and right tailed = greater
#' @param alternative Specify alternative, left tailed = "less" and right tailed = "greater"
#' @param n.burnin Number of burn-in iterations in the MCMC.
#' @param n.iter Total number of iterations in the MCMC. Post burn-in = n.iter-n.burnin
#' @return summary from JAGS output containing the estimates of all parameters in the Bayesian hierarchical model.
#' @export
#' @importFrom R2jags jags
#' @importFrom assertthat assert_that
ppos_surrogate2 <- function(Hy,
                     Hy.se,
                     rho.hat,
                     tau.dist,
                     tau.prior,
                     alpha,
                     X.cs,
                     sigma.fs,
                     sigma.cs,
                     alternative ,
                     n.burnin=4000, n.iter=10000
){


  assert_that(is.matrix(Hy))
  assert_that(ncol(Hy)==2)

  assert_that(is.matrix(Hy.se))
  assert_that(ncol(Hy.se)==2)

  data <- data.prep(Hy,Hy.se,rho.hat,alpha,sigma.fs,sigma.cs,alternative,X.cs,tau.prior)
  para.save=c("eta","ppos","Sigma","theta","rho","Tf")

  prior <- paste(tau.dist)
  model.file0 <-  system.file("model", paste0(prior, ".txt"), package = "ppos")

  #print(getwd())
  jags.mod.fit <- jags(data = data, inits = NULL,
                       parameters.to.save = para.save, n.chains = 2, n.iter = n.iter,n.thin = 1,
                       n.burnin = n.burnin, model.file = model.file0)
  Rhat.max <- max(jags.mod.fit[["BUGSoutput"]][["summary"]][,"Rhat"])
  if(Rhat.max > 1.1)
    warning("Maximal Rhat > 1.1. Consider increasing n.burnin and n.iter MCMC parameter.")
  out <- as.data.frame(jags.mod.fit$BUGSoutput$summary)
  return(out)

}


#' Posterior analysis: using MAP prior
#'
#' @param Hy Data set containing treatment effects from historical studies.
#' @param Hy.se Data set (historical) containing square of standard error of treatment effects from historical studies.
#' @param tau.dist default = NULL
#' @param tau.prior Parameters of prior distribution for \code{tau}. default=NULL estimates \code{tau} using DerSimonian and Laird (DL) approach.
#' @param alpha Two-sided alpha.
#' @param sigma.fs Estimated variance of the treatment effect yet to be observed.
#' @param alternative Specify alternative, left tailed = "less" and right tailed = "greater"
#' @param n.burnin Number of burn-in iterations in the MCMC.
#' @param n.iter Total number of iterations in the MCMC. Post burn-in = n.iter-n.burnin
#' @return summary from JAGS output containing the estimates of all parameters in the Bayesian hierarchical model.
#' @export
ppos_MAP <- function(Hy,
                      Hy.se,
                      tau.dist,
                      tau.prior,
                      alpha,
                      sigma.fs,
                      alternative ,
                      n.burnin=4000, n.iter=10000
){

  if(is.null(tau.dist) & is.null(tau.prior)){
    print("Estimating tau2 using DerSimonian and Laird (DL) approach")
    data.map <- data.prep.mapDL(Hy,Hy.se)
    data.map$alternative <- ifelse(alternative=="less",1,0)
  }
  #path <- path.package('ppos')
  #file <- file.path(path, 'model', 'MAPmodel2_ppos.txt')
  model.file <-  system.file("model","MAPmodel2_ppos.txt", package = "ppos")
  jags.mod.fit <- MAPrun(model.file0=model.file,data=data.map,sigma.fs,alpha,alternative,n.burnin=4000, n.iter=10000)
  Rhat.max <- max(jags.mod.fit[["BUGSoutput"]][["summary"]][,"Rhat"])
  if(Rhat.max > 1.1)
    warning("Maximal Rhat > 1.1. Consider increasing n.burnin and n.iter MCMC parameter.")
  out <- as.data.frame(jags.mod.fit$BUGSoutput$summary)
  return(out)
}



