#' Decision Boundary to obtain the critical region R
#'
#' Calculates the decision boundary for a 1 sample design. This is the
#' critical value at which the decision function will change from 0
#' (failure) to 1 (success).
#'
#' @param alpha Two sided type I error.
#' @param sigma Estimated sigma for the future outcome \eqn{T_f}. In the case of log hazard ratio it can be considered as \eqn{\sqrt(4/d)}, where \eqn{d} are the number of events.
#' @param alternative Specify alternative (left tailed or right tailed)
#' @returns critical value
#' @importFrom stats qnorm
decision_boundary <- function(alpha,
                              sigma,
                              alternative = c("less", "greater")){
  if(alternative=="less"){
    Za <- qnorm(alpha/2)
  }else if(alternative=="greater"){
    Za <- qnorm(1-(alpha/2))
  }

  delta <- sigma *Za
  return(delta)
}




#' This function converts the within-study standard errors obtained from historical studies to variance-covariance or precision matrices required as an input
#' for jags model fit given the within-study correlartion coefficint estimate
#'
#' @param sigi a \eqn{K \times 1} within-study standard errors corresponding to primary outcome
#' @param deltai a \eqn{K \times 1} within-study standard errors corresponding to surrogate/secondary outcome
#' @param rho.hat user input of within-study correlation coefficient estimate
#' @returns a 3-d array of dimension \eqn{K \times m \times m} of precision matrices
#' @returns a 3-d array of dimension \eqn{m \times m \times K} of variance-covaraince matrices

input.data.prec.bhm <- function(sigi,deltai,rho.hat){

  rho=rho.hat
  TT <- array(dim=c(2,2,length(deltai)))
  for(i in 1:length(deltai)){
    cov <- sigi[i]*deltai[i] * rho
    TT[,,i] <- matrix(c(sigi[i]^2, cov,cov, deltai[i]^2  ),ncol=2,byrow = T)
  }

  mat.tt <- lapply(1:dim(TT)[3],function(ii) solve(TT[,,ii]))


  test.arr <- array(dim=c(2,2,length(mat.tt)))
  for(i in 1:length(mat.tt)){
    test.arr[,,i] <- mat.tt[[i]]
  }


  # Convert R-3D array to BUGS 3-D array
  dim <- dim(test.arr)
  test.arr.winB <- array(0,dim=c(dim[3],dim[1],dim[2]))

  for(i in 1:dim[1]){
    for(j in 1:dim[2]){
      for(k in 1:dim[3]){
        test.arr.winB[k,i,j] <-  test.arr[i,j,k]
      }
    }
  }
  return(list(prec.mat=test.arr.winB,var.mat=TT))
}


data.prep <- function(Hy,
                      Hy.se,
                      rho.hat,
                      alpha,
                      sigma.fs,
                      sigma.cs,
                      alternative,
                      X.cs,
                      tau.prior)
  {
  # update dy and dy.se for current trial
  dy <- rbind(Hy,c(NA,X.cs))
  dy.se <-  rbind(Hy.se,c(sigma.fs,sigma.cs))

  prec.hat <- input.data.prec.bhm(sigi=dy.se[,1],deltai=dy.se[,2],rho.hat=rho.hat)

  delta <- decision_boundary(alpha,
                             sigma.fs,
                        alternative)
  n.r <- nrow(dy)
  dat = list(
    I=n.r,
    thet.hat = dy,
    prec.hat = prec.hat$prec.mat[1:(n.r-1),,],
    delta = delta,
    Sigma.curr = solve(prec.hat$prec.mat[(n.r),,]),
    tau.prior=tau.prior,
    alternative =ifelse(alternative=="less",1,0)
  )
  return(dat)
}

se_from_CI <- function(upper.limit,lower.limit,L){
  se=(upper.limit-lower.limit)/(2*qnorm(1-(1-L/100)/2))
  return(se)
}


## ================================
#  MAP model with ppos
## ================================
#' @importFrom meta metamean
data.prep.mapDL <- function(Hy,Hy.se){
  data <- NULL
  data$X.k <- Hy
  data$sigma2.k <- Hy.se^2
  # Hyper parameters for the overall mean
  # theta ~ N(t0,sig2.0), hyperprior on theta=overall mean
  t0=0
  sig2.0=10^6
  data.MAP <- c(data,t0=t0,tau2.0=1/sig2.0,K=length(data$X.k))

  M1 <- metamean(
    n=rep(1,data.MAP$K),
    mean=data.MAP$X.k,
    sd=sqrt(data.MAP$sigma2.k),
    method.tau= ("DL"))
  tau2.DS = 1/M1$tau2
  data.MAP$tau2.DS <- tau2.DS
  return(data.MAP)
}


MAPrun <- function(model.file0,data,sigma.fs,alpha,alternative,n.burnin, n.iter){
  data.MAP <- data

  #
  para.save = c("ppos","t","Tf","t.star","X.k","CPOinv")

  model.data2 = data.MAP
  model.data2$K <- model.data2$K + 1
  model.data2$X.k <- c(model.data2$X.k,NA)
  model.data2$sigma2.k <- c(model.data2$sigma2.k,sig2.f=sigma.fs^2)

  delta <- decision_boundary(alpha,
                             sigma.fs,
                             alternative)

  model.data2$delta <- delta #qnorm(alph)*sqrt(sig2.f)

  MAP.mod.fit2 <- jags(data = model.data2, inits = NULL,
                       parameters.to.save = para.save, n.chains = 2, n.iter = n.iter,n.thin = 1,
                       n.burnin = n.burnin, model.file = model.file0)


  # CPO1 <- 1/MAP.mod.fit2$BUGSoutput$mean$CPOinv
  # LPML1 <- sum(log(CPO1))
  # LPML1
  # MAP.mod.fit2$BUGSoutput$DIC
  # #print(MAP.mod.fit2)
  # tf.map <- MAP.mod.fit2[["BUGSoutput"]][["sims.list"]][["Tf"]]
  # summary.model2 <- MAP.mod.fit2[["BUGSoutput"]][["summary"]]
  # summary.model2<- data.frame(summary.model2)
  #
  #
  # theta.star.mcmc02 <- MAP.mod.fit2$BUGSoutput$sims.list$t.star[,-model.data2$K] #[c(n.burnin:n.iter),]
  # fx.mcmc <- sapply(1:nrow(theta.star.mcmc02), function(ii) fx_iter(data,ii,theta.star.mcmc02))
  # fx.mcmc <- t(fx.mcmc)
  # dim(fx.mcmc)
  #
  # dd02 <- data.frame(model="MAP-DL",M=1,deviance=MAP.mod.fit2$BUGSoutput$summary["deviance",1],
  #                    LPML=LPML(fx.mcmc),LPML1=LPML1,
  #                    DIC=MAP.mod.fit2$BUGSoutput$DIC,
  #                    tau2=1/data.MAP$tau2.DS,
  #                    ppos=summary.model2["ppos",1],
  #                    Tf=summary.model2["Tf",c("mean","sd","X2.5.","X50.","X97.5.")])
  # colnames(dd02)[11:13] <- c("Tf.2.5." ,"Tf.50.", "Tf.97.5.")
  # model.perf.MAP2 <- dd02
  # OUT<- NULL
  # OUT$model.perf <- model.perf.MAP2
  # OUT$summary.model <- summary.model2
  return(MAP.mod.fit2)
}





#' @importFrom dplyr %>%
#' @importFrom dplyr filter
#' @importFrom stringr str_detect
#' @importFrom forestplot forestplot
#' @importFrom forestplot fp_set_style

# THIS PLOT IS FOR BIVARAITE CASE ----
forestplot2 <- function(out,primary,surr){
  summary.out <- cbind(out, names=rownames(out))
  ind.theta <- summary.out %>% filter(str_detect(names, "^theta"))
  thet1 <- ind.theta %>% filter(str_detect(names, "1]$"))
  thet2 <- ind.theta %>% filter(str_detect(names, "2]$"))
  eta1 <- summary.out["eta[1]",]
  eta2 <- summary.out["eta[2]",]
  thet1 <-data.frame(rbind(thet1,eta1))
  thet2 <-data.frame(rbind(thet2,eta2))
  # forest plot
  plot.delta1 = data.frame(mean = thet1$mean, lower = thet1$X2.5., upper = thet1$X97.5.)
  tabletext = c(paste("Study",1:(nrow(plot.delta1)-1)),"MAP mean")
  plot.delta2 = data.frame(mean = thet2$mean, lower = thet2$X2.5., upper = thet2$X97.5.)
  if(surr=="PFS"){
    tit2 <- paste("log(",surr,")",sep="")
  }
  if(surr=="ORR"){
    tit2 <- paste(expression(Delta),surr)
  }
  if(primary %in% c("PFS","OS")){
    tit1 <- paste("log(",primary,")",sep="")
  }

  fp1 <- forestplot(tabletext,
                    plot.delta1,title=tit1)
  fp2 <- forestplot(tabletext,
                    plot.delta2,title=tit2)

  return(list(fp1=fp1,fp2=fp2))
}

# THIS PLOT IS FOR UNIVARIATE CASE ----
#' Forest plot of the historical data. The K+1th study shows the predicted value of the future tretament effect.
#'
#' @param out Output from JAGS from function "ppos_MAP"
#' @return forest plot
#' @export
forestplot3 <- function(out){
  summary.out <- cbind(out, names=rownames(out))
  ind.theta <- summary.out %>% filter(str_detect(names, "^t.star"))
  eta <- summary.out["t",]
  thet1 <-data.frame(rbind(ind.theta,eta))

  # forest plot
  plot.delta1 = data.frame(mean = thet1$mean, lower = thet1$X2.5., upper = thet1$X97.5.)
  tabletext = c(paste("Study",1:(nrow(plot.delta1)-1)),"MAP mean")
  clip1 <- min(unlist(plot.delta1)) - 0.05
  clip2 <- max(unlist(plot.delta1)) + 0.05
  ticks <- seq(-1,1,by=0.2)
  fp1 <- forestplot(tabletext,
                    mean=plot.delta1$mean,
                    lower=plot.delta1$lower,
                    upper=plot.delta1$upper,
                    clip=c(clip1,clip2),
                    xticks = ticks,#c(clip1, -0.5, 0,0.5,clip1),
                    title=" ") |>
    fp_set_style(box = "royalblue",
                 lines = "darkblue",
                 summary = "royalblue")

  return(list(fp1=fp1))
}

#' Predictive probability of success PPoS
#'
#' @param jags.output Output from JAGS
#' @return ppos value
#' @export
ppos <- function(jags.output){
  res <- round(jags.output["ppos",1],3)
  return(res)
}
