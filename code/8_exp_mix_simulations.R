#############################################################################################
# Variational asymptotics and sandwich estimation
# 
# File: exp_mix_sim_replication.R
# -----------------------------
# Replicates simulation study and makes plots for the mixture of exponentials model
#
#############################################################################################

source('code/exp_mixture_vb_fns.R')

# Define simulation parameters
b <- 1000
cases <- expand.grid(theta=1,
                     lambda=2,
                     n=1000,
                     Xdist=c("exp", "gamma"),
                     Zdist="exp")

# Function to run one simulation
sim_function <- function(case) {
  print(case)
  n <- case$n
  # Generate all data for this case
  if(case$Xdist == "exp") X <- rexp(n*b, rate=case$theta)
  if(case$Xdist == "gamma") X <- rgamma(n*b, shape=3, rate=3*case$theta)
  Z <- rexp(n*b, rate=case$lambda)
  Y <- rexp(n*b, rate=Z * case$theta)
  Xmat <- matrix(X, nrow=b)
  Zmat <- matrix(Z, nrow=b)
  Ymat <- matrix(Y, nrow=b)
  
  theta0 <- case$theta
  lambda0 <- case$lambda
  
  require(plyr)
  ldply(1:b, function(bi) {
    if(bi %% 100 == 0) print(bi)
    X <- Xmat[bi,]
    Y <- Ymat[bi,]
    Z <- Zmat[bi,]
    N <- n
    
    # Initialize parameters for VB algorithm
    # Note: the intialization is very important. VB has many local modes, especially as the 
    #   sample size increase. Empirically it appears that all of these modes are consistent,
    #   but may have different relative efficiencies. For instance, initializing arbitrarily 
    #   leads to a non-efficient mode, and information/sandwich-based CIs will be a little
    #   too small, even asymptotically. Initializing at the truth leads to an unrealistically
    #   good estimate, and CIs will be conservative. Here we use a simple moment-based initialization
    #   that appears to get the algorithm to the true maximizer of the ELBO.
    theta_start <- 1/ mean(X)
    lambda_start <- theta_start * median(Y) 
    Z_start <- 2 / (theta_start * Y + lambda_start)
    
    # Estimate mean-field VB
    est <- exp_mixture_vb(X, Y, Z_start=Z_start, theta_start = theta_start, lambda_start = lambda_start)
    thetahat <- est$alphatheta/est$betatheta
    lambdahat <- est$alphalambda/est$betalambda
    
    # Get information & sandwich covariances
    covs <- exp_mixture_cov(X, Y, thetahat =thetahat, lambdahat = lambdahat, etaZ=est$alphaZ/est$betaZ, alphaZ = est$alphaZ)
    
    # Construct CIs
    vb_theta <- qgamma(c(.025, .975), shape=est$alphatheta, rate=est$betatheta)
    vb_lambda <- qgamma(c(.025, .975), shape=est$alphalambda, rate=est$betalambda)
    vb_theta_joint <- qgamma(c((1-sqrt(.95))/2,(1+sqrt(.95))/2),shape=est$alphatheta, rate=est$betatheta)
    vb_lambda_joint <- qgamma(c((1-sqrt(.95))/2,(1+sqrt(.95))/2),shape=est$alphalambda, rate=est$betalambda)
    
    sand_theta <- qnorm(c(.025, .975), thetahat, sqrt(covs$sand_cov[1,1]/n)) 
    sand_lambda <- qnorm(c(.025, .975), lambdahat, sqrt(covs$sand_cov[2,2]/n))
    sand_joint_thetalambda <- as.vector(t(c(theta0, lambda0) - c(thetahat, lambdahat)) %*% solve(covs$sand_cov/n) %*% (c(theta0, lambda0) - c(thetahat, lambdahat))) <= qchisq(.95, 2)
    
    df <- data.frame(n=n, theta=case$theta, lambda=case$lambda, Xdist=case$Xdist, Zdist=case$Zdist,
                     vb_contains_theta=vb_theta[1] <= theta0 & theta0 <= vb_theta[2],
                     vb_contains_lambda=vb_lambda[1] <= lambda0 & lambda0 <= vb_lambda[2],
                     vb_contains_joint=vb_theta_joint[1] <= theta0 & theta0 <= vb_theta_joint[2] & vb_lambda_joint[1] <= lambda0 & lambda0 <= vb_lambda_joint[2],
                     sand_contains_theta=sand_theta[1] <= theta0 & theta0 <= sand_theta[2],
                     sand_contains_lambda=sand_lambda[1] <= lambda0 & lambda0 <= sand_lambda[2],
                     sand_contains_joint=sand_joint_thetalambda,
                     vb_theta_var=est$alphatheta / est$betatheta^2,
                     vb_lambda_var=est$alphalambda / est$betalambda^2,
                     vb_joint_cov=0,
                     sand_theta_var=covs$sand_cov[1,1]/n,
                     sand_lambda_var=covs$sand_cov[2,2]/n,
                     sand_joint_cov=covs$sand_cov[1,2]/n,
                     thetahat=thetahat,
                     lambdahat=lambdahat)
    return(df)
  })
}

# Run simulations
set.seed(442384)
library(plyr)
results <- adply(cases, 1, sim_function, .progress="text")

save(results, file='data/exp_mix_sim_results.Rdata')

load('exp_mix_sim_results.Rdata')
name_parts <- strsplit(names(results), "_")

# FREQUENTIST COVERAGE
# Get relevant columns 
cont <- which(unlist(lapply(name_parts, function(l) "contains" %in% l)))

# Compute coverages
coverages <- ddply(results[,1:14], .(n,theta, lambda, Xdist, Zdist), function(df) {
  cis <- t(sapply(colSums(df[,cont]), function(col) c(binom.test(col, n=nrow(df))$conf.int)))
  data.frame(Covariance=unlist(lapply(name_parts[cont], function(l) l[1])), Parameter=unlist(lapply(name_parts[cont], function(l) l[3])), Coverage = colMeans(df[,cont]), CoverageLower=cis[,1], CoverageUpper=cis[,2])
})

coverages$Parameter <- factor(coverages$Parameter, levels=c("theta", "lambda", "joint"), labels=c("theta", "lambda", "(list(theta, lambda))"))
coverages$Covariance <- factor(coverages$Covariance, levels=c("vb", "info", "sand"), labels=c("VB", "Info.", "Sand."))

library(xtable)
xtable(coverages[,c(6,7,8)])
