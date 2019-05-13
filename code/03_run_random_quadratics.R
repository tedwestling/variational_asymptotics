###############################################################################
# REPLICATION FILES FOR:
# Beyond prediction: A framework for inference with variational approximations in mixture models
# Ted Westling and Tyler McCormick
#
# file: 03_run_random_quadratics.R
# 
# The file 01_compile_data.R needs to have already been run, which will have created
# the file data/random_quadratic_data.Rdata,
#
# This script computes the various estimates of the logistic random quadratic model
# for the marijuana use variables of the NLSY97 data. There are six fixed
# effects corresponding to the intercept, slope, and quadratic effect of time
# for the average male and female separately.
# 
# Author: Ted Westling (ted.westling@gmail.com)
###############################################################################

## NOTE : This code was developed using version 1.1-7 of lme4
#         Newer versions of lme4 have changed certain package syntax and warnings,
#         and running the code with the current version of lme4 may produce warnings or errors
#  To install version 1.1-7 of lme4, run the following code, then restart R.
# library(devtools)
# install_version("lme4", version = "1.1-7", repos = "http://cran.us.r-project.org")

library(MASS)
library(numDeriv)
library(plyr)
library(fastGHQuad)
library(lme4)
library(mvQuad)
library(parallel)

rule100 <- gaussHermiteData(100)

source('code/mixed_logit_functions.R')
nodes <- 1 # we ran this file on a personal computer with no parallelization

# Load compiled data
# The file 01_compile_data.R needs to have already been run, which will have created
# the file data/random_quadratic_data.Rdata,
load('data/random_quadratic_data.Rdata')

# For code testing, uncomment the following line to take a subsample
# data <- data[1:1000]


beta_length <- ncol(data[[1]]$x)
rand_eff_length <- 3
Nsub <- length(data)

# Run glmer estimate
x_mat <- do.call("rbind", lapply(data, function(d) d$x))
y_mat <- unlist(lapply(data, function(d) d$y))
z_mat <- do.call("rbind", lapply(data, function(d) d$z))
sub_mat <- rep(1:Nsub, times=unlist(lapply(data, function(d) nrow(d$x))))

# Note: The following command may produce a warning regarding non-convergence of the optimization used to estimate
#       the GLMM. We are only using this estimate to get a starting point for the variational algorithm below, 
#       not to report.
glmer_est <- glmer(y_mat ~ x_mat[,2] + x_mat[,3] + x_mat[,4] + x_mat[,5] + x_mat[,6] + (1 + z_mat[,2] + z_mat[,3] | sub_mat), family='binomial', verbose=1, control=glmerControl(check.conv.grad="ignore", check.conv.singular="ignore", check.conv.hess="ignore"))

glmer_beta <- summary(glmer_est)$coefficients[,1]

# The way we defined the model was in terms of separate intercepts for males and females, with no global intercept,
# so we have to  add the intercept for the coefficient for males
glmer_beta[4] <- glmer_beta[1] + glmer_beta[4]

glmer_cor <- attr(summary(glmer_est)$varcor[[1]], "correlation")
glmer_sd <- attr(summary(glmer_est)$varcor[[1]], "stddev")
glmer_Sigma <- diag(glmer_sd) %*% glmer_cor %*% diag(glmer_sd)

glmer_Sigma_D <- log(LDL(glmer_Sigma)$D)
glmer_Sigma_L <- LDL(glmer_Sigma)$L
glmer_theta <- c(glmer_beta, glmer_Sigma, glmer_Sigma_D, glmer_Sigma_L[lower.tri(glmer_Sigma_L, diag=FALSE)])

# Run variational EM algorithim
set.seed(187112)
beta_start <- glmer_beta + rnorm(length(glmer_beta))
Sigma_start <- diag(rep(1,3))

# Note: the following line of code takes approximately 12 hours to complete because we set the
# tolerance very low and there are over 100K observations -- we ran it overnight
vem_opt <- mixed_logit_vem(data=data, beta_start=glmer_beta, Sigma_start=diag(rep(1,3)), verbose=TRUE, nodes=nodes, tol=.0001)

opt_par <- par_vec_to_list(c(vem_opt$beta, unlist(vem_opt$psi)), K=beta_length, P=3, n=length(data))

vem_beta <- opt_par$beta
vem_Sigma <- opt_par$Sigma
vem_Sigma_D <- log(LDL(vem_Sigma)$D)
vem_Sigma_L <- LDL(vem_Sigma)$L
vem_theta <- c(vem_beta, vem_Sigma_D, vem_Sigma_L[lower.tri(vem_Sigma_L, diag=FALSE)])

# Compute sandwich covariance
# Note: the following line of code takes several hours to complete
sand_cov <- sandwich_cov(beta=vem_beta, Sigma_D=vem_Sigma_D, Sigma_L=vem_Sigma_L, m=opt_par$m, D=opt_par$D, L=opt_par$L, data=data, verbose=TRUE, nodes=nodes)

# Compute score and information
grid <- createNIGrid(dim=3, type='GHN', level=10)
# The next command takes about 10 minutes to run
vem_score_info <- score_info(vem_beta, vem_Sigma_D, vem_Sigma_L, data=data, verbose=TRUE, nodes=nodes)

os_theta <- vem_theta + solve(vem_score_info$obs_info) %*% vem_score_info$score
os_beta <- os_theta[1:6]
os_Sigma_D <- os_theta[7:9]
os_Sigma_L <- diag(rep(1,3))
os_Sigma_L[lower.tri(os_Sigma_L, diag=FALSE)] <- os_theta[10:12]
os_Sigma <- os_Sigma_L %*% diag(exp(os_Sigma_D)) %*% t(os_Sigma_L)

# Again, about 10 minutes
os_score_info <- score_info(os_beta, os_Sigma_D, os_Sigma_L, data=data, verbose=TRUE, nodes=nodes)

n_sim <- 100000
# For code testing, set n_sim to 1000
# n_sim <- 1000
set.seed(99099092)
sub_effs <- mvrnorm(n_sim, rep(0,3), vem_Sigma)

# Note: the following code takes several hours to run
start_time <- Sys.time()
prof_func_hats <- mclapply(1:n_sim, function(i) {
  if(i %% 500 == 0) {
    time_left <- (c(difftime(Sys.time(), start_time, units="secs")) / i) * (n_sim - i)
    cat("Obs", i, "of", n_sim, ",", round(c(time_left)/60, 2), "minutes left\n")
  }
  samp <- sample(1:length(data), 1)
  z <- data[[samp]]$z
  x <- data[[samp]]$x
  prob <-  c(1/(1+exp(-(x %*% vem_beta + z %*% sub_effs[i,]))))
  y <- rbinom(length(prob), 1, prob)
  this_opt <- optim(rep(0, 9), fn=one_elbo_theta0, gr=one_elbo_grad_theta0, beta=vem_beta, Sigma=vem_Sigma, Yi=y, Zi=z, Xi=x, control=list(fnscale=-1), method='L-BFGS-B')
  mi <- this_opt$par[1:3]
  Di <- this_opt$par[4:6]
  Li <- diag(rep(1,3))
  Li[lower.tri(Li, diag=FALSE)] <- this_opt$par[7:9]
  grad <- one_elbo_thetagrad(beta=vem_beta, Sigma=vem_Sigma, mi=mi, Di=Di, Li=Li, Yi=y, Zi=z, Xi=x)
  return(list(this_opt, grad))
}, mc.cores=nodes)

save(glmer_est, vem_opt, sand_cov, vem_score_info, os_score_info, prof_func_hats, file='estimates/rand_quadratics_estimate.Rdata')

