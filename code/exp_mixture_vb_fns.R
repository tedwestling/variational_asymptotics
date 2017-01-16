ensure_package <- function(package_name) {
  for (i in 1:length(package_name)) {
    name <- package_name[i]
    if(!require(name, character.only=TRUE, quietly=TRUE)) {
      install.packages(name, repos='http://cran.cnr.Berkeley.edu')
      library(name, character.only=TRUE)
    }
  }
}

ensure_package('ggplot2')
ensure_package('plyr')
ensure_package('MASS')
ensure_package('mgcv')
ensure_package('reshape')
ensure_package('numDeriv')
ensure_package('boot')
ensure_package('grid')
ensure_package('gridExtra')
ensure_package('Matrix')

exp_mixture_cov <- function(X, Y, thetahat, lambdahat, etaZ, alphaZ) {
  dt <- 2/thetahat - X - etaZ * Y
  dtt <- rep(-2/thetahat^2, length(X))
  dl <- 1/lambdahat - etaZ
  dll <- rep(-1/lambdahat^2, length(X))
  dn <- -thetahat * Y - lambdahat + 2/etaZ
  dnn <- -2/etaZ^2
  dtn <- -Y
  dln <- rep(-1, length(Y))
  dtl <- 0
  
  dndt <- -Y * etaZ^2 /2
  dndl <- -etaZ^2 / 2
  
  grads <- cbind(dt, dl)
  grad_xprods <- cbind(grads[,1]^2, grads[,1] * grads[,2], grads[,1] * grads[,2], grads[,2]^2)
  Bhat <- matrix(colMeans(grad_xprods), nrow=2)
  
  hesss <- matrix(colMeans(cbind(dtt + dtn * dndt, dtl + dln * dndt, dtl + dln * dndt, dll + dln * dndl)), nrow=2)
  
  Ahat_inv <- solve(hesss)
  Bhat <- matrix(colMeans(grad_xprods), nrow=2)
  info_cov <- -matrix(c( -2 * lambdahat * mean((2*thetahat*Y + lambdahat)/(Y * thetahat^2 + lambdahat*thetahat)^2), 2 * mean(Y / (Y * thetahat + lambdahat)^2),
                              2 * mean(Y / (Y * thetahat + lambdahat)^2), -lambdahat^{-2} + (2 / thetahat^2) * mean(1/(Y + lambdahat / thetahat)^2)), nrow=2, ncol=2)
  list(sand_cov = Ahat_inv %*% Bhat %*% Ahat_inv, inv_info_cov = solve(info_cov))
}

exp_mixture_vb <- function(X, Y, alpha0theta=1, beta0theta=1, alpha0lambda=1, beta0lambda=1, theta_start=NULL, lambda_start=NULL, Z_start=NULL, verbose=FALSE) {
  xsum <- sum(X)
  n <- length(X)
  alphatheta <- 2*n + alpha0theta
  alphalambda <- n + alpha0lambda
  alphaZ <- 2
  if(is.null(theta_start)) betatheta <- alphatheta * mean(X)
  else betatheta <- alphatheta / theta_start
  if(is.null(Z_start)) betaZ <- alphaZ * alphatheta * Y / betatheta
  else betaZ <- alphaZ / Z_start
  if(is.null(lambda_start)) betalambda <- sum(alphaZ/betaZ) + beta0lambda
  else betalambda <- alphalambda / lambda_start
  elbo <- (2*n + alpha0theta - 1) * (digamma(alphatheta) - log(betatheta)) - (beta0theta + xsum)* alphatheta/betatheta  +
    (n + alpha0lambda - 1) * (digamma(alphalambda) - log(betalambda)) - (beta0lambda + sum(alphaZ/betaZ)) * alphalambda / betalambda +
    sum(digamma(alphaZ) - log(betaZ)) - (alphatheta / betatheta) * sum(alphaZ * Y / betaZ) -
    gamma_entropy(alphatheta, betatheta) - gamma_entropy(alphalambda, betalambda) - sum(gamma_entropy(alphaZ, betaZ))
  elbo_last <- NULL
  
  while(is.null(elbo_last) || abs((elbo_last - elbo)/elbo_last) > .001) {
    if(verbose) {
      print(paste("theta =", alphatheta/betatheta, "; lambda =", alphalambda/betalambda, "; elbo =", round(elbo, 4)))
    }
    betatheta <- xsum + sum(alphaZ*Y/betaZ) + beta0theta
    betalambda <- sum(alphaZ/betaZ) + beta0lambda
    betaZ <- Y * alphatheta / betatheta + alphalambda/betalambda
    elbo_last <- elbo
    elbo <- (2*n + alpha0theta - 1) * (digamma(alphatheta) - log(betatheta)) - (beta0theta + xsum)* alphatheta/betatheta  +
      (n + alpha0lambda - 1) * (digamma(alphalambda) - log(betalambda)) - (beta0lambda + sum(alphaZ/betaZ)) * alphalambda / betalambda +
      sum(digamma(alphaZ) - log(betaZ)) - (alphatheta / betatheta) * sum(alphaZ * Y / betaZ) -
      gamma_entropy(alphatheta, betatheta) - gamma_entropy(alphalambda, betalambda) - sum(gamma_entropy(alphaZ, betaZ))
  }
  return(list(alphatheta=alphatheta, betatheta=betatheta, alphalambda=alphalambda, betalambda=betalambda, alphaZ=alphaZ, betaZ=betaZ, elbo=elbo))
}

gamma_entropy <- function(alpha, beta) alpha - log(beta) + lgamma(alpha) + (1-alpha) * digamma(alpha)

exp_mixture_gibbs <- function(X, Y, alpha0theta=1, beta0theta=1, alpha0lambda=1, beta0lambda=1, theta_start=NULL, lambda_start=NULL, Z_start=NULL, verbose=FALSE, burnin=1000, samples=1000, thin=10) {
  xsum <- sum(X)
  n <- length(X)
  if(is.null(theta_start)) theta_start <- 1/ mean(X)
  if(is.null(lambda_start)) lambda_start <- theta_start * median(Y) 
  if(is.null(Z_start)) Z_start <- 2 / (theta_start * Y + lambda_start)
  
  tot_samp <- burnin + samples * thin
  par_t <- list(Z=Z_start, theta=theta_start, lambda=lambda_start)
  samples <- list(theta=rep(NA, tot_samp), lambda=rep(NA, tot_samp), Z=matrix(NA, nrow=tot_samp, ncol=n))
  for(t in 1:tot_samp) {
    if(verbose) print(t)
    samples$theta[t] <- par_t$theta <- rgamma(1, alpha0theta + 2*n, beta0theta + xsum + sum(Y * par_t$Z))
    samples$lambda[t] <- par_t$lambda <- rgamma(1, alpha0lambda + n, beta0lambda + sum(par_t$Z))
    for(i in 1:n) samples$Z[t,i] <- par_t$Z[i] <- rgamma(1, 2, par_t$theta * Y[i] + par_t$lambda)
  }
  return(samples)
}


exp_mixture_ln_vb <- function(X, Y, alpha0theta=1, beta0theta=1, alpha0lambda=1, beta0lambda=1, theta_start=NULL, lambda_start=NULL, Z_start=NULL, verbose=FALSE) {
  xsum <- sum(X)
  n <- length(X)
  alphatheta <- 2*n + alpha0theta
  alphalambda <- n + alpha0lambda
  sigmaZ <- sqrt(1/2)
  if(is.null(theta_start)) betatheta <- alphatheta * mean(X)
  else betatheta <- alphatheta / theta_start
  if(is.null(Z_start)) muZ <- log(betatheta/(alphatheta * Y)) - sigmaZ^2 / 2
  else muZ <- log(Z_start) - sigmaZ^2 / 2
  if(is.null(lambda_start)) betalambda <- sum(exp(muZ + sigmaZ^2/2)) + beta0lambda
  else betalambda <- alphalambda / lambda_start
  elbo <- (2*n + alpha0theta - 1) * (digamma(alphatheta) - log(betatheta)) - (beta0theta + xsum)* alphatheta/betatheta  +
    (n + alpha0lambda - 1) * (digamma(alphalambda) - log(betalambda)) - (beta0lambda + sum(exp(muZ + sigmaZ^2/2))) * alphalambda / betalambda +
    sum(muZ) - (alphatheta / betatheta) * sum( Y * exp(muZ + sigmaZ^2/2)) -
    gamma_entropy(alphatheta, betatheta) - gamma_entropy(alphalambda, betalambda) - sum(ln_entropy(muZ, sigmaZ))
  elbo_last <- NULL
  
  while(is.null(elbo_last) || abs((elbo_last - elbo)/elbo_last) > .001) {
    if(verbose)  print(paste("theta =", alphatheta/betatheta, "; lambda =", alphalambda/betalambda, "; elbo =", round(elbo, 4)))
    betatheta <- xsum + sum(Y * exp(muZ + sigmaZ^2/2)) + beta0theta
    betalambda <- sum(exp(muZ + sigmaZ^2/2)) + beta0lambda
    muZ <- - sigmaZ^2 / 2 + log(2 / (Y * alphatheta / betatheta + alphalambda/betalambda))
    elbo_last <- elbo
    elbo <- (2*n + alpha0theta - 1) * (digamma(alphatheta) - log(betatheta)) - (beta0theta + xsum)* alphatheta/betatheta  +
      (n + alpha0lambda - 1) * (digamma(alphalambda) - log(betalambda)) - (beta0lambda + sum(exp(muZ + sigmaZ^2/2))) * alphalambda / betalambda +
      sum(muZ) - (alphatheta / betatheta) * sum( Y * exp(muZ + sigmaZ^2/2)) -
      gamma_entropy(alphatheta, betatheta) - gamma_entropy(alphalambda, betalambda) - sum(ln_entropy(muZ, sigmaZ))
  }
  return(list(alphatheta=alphatheta, betatheta=betatheta, alphalambda=alphalambda, betalambda=betalambda, muZ=muZ, sigmaZ=sigmaZ, elbo=elbo))
}

ln_entropy <- function(mu, sigma) log(sigma) + mu
