###############################################################################
# REPLICATION FILES FOR:
# Beyond prediction: A framework for inference with variational approximations in mixture models
# Ted Westling and Tyler McCormick
#
# file: 09_timing_simulations.R
#
# The file 01_compile_data.R needs to have already been run, which will have created
# the file data/random_quadratic_data.Rdata.
#
# This script runs a single timing simulation using the logistic regression 
# with random quadratics simulation.
# The true parameters were drawn from the estimates for the NLSY97 data.
# MLE (with glmer), variational EM, variational EM + one step, and MLE are used
# to estimate the parameters. Here is the slurm code used to run all 1000 
# simulations:
# !/bin/bash
# #SBATCH --job-name timing_sim
# #SBATCH --partition medium     # Slurm partition to use
# #SBATCH --ntasks 1         # Number of tasks to run. By default, one CPU core will be allocated per task
# #SBATCH --cpus-per-task=1
# #SBATCH --time 2-12:00        # Wall time limit in D-HH:MM
# #SBATCH --mem-per-cpu=800     # Memory limit for each tasks (in MB)
# #SBATCH -o outs/timing_sim_%A_%a.out    # File to which STDOUT will be written
# #SBATCH -e outs/timing_sim_%A_%a.err    # File to which STDERR will be written
# #SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
# #SBATCH --mail-user=tgwest@uw.edu # Email to which notifications will be sent
# 
# R CMD BATCH --vanilla --no-restore --no-save --args n ${1} < /code/9_timing_simulations.R  > outs/timing_sim_n${1}_"${SLURM_ARRAY_TASK_ID}".Rout
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
library(batch)
library(lme4)
library(mvQuad)
library(parallel)
library(DescTools)

rule100 <- gaussHermiteData(100)

source('code/mixed_logit_functions.R')
nodes <- 1
task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

parseCommandArgs()


# Load compiled data
# The file 1_compile_data.R needs to have already been run, which will have created
# the file data/random_quadratic_data.Rdata.
load('data/random_quadratic_data.Rdata')

true_beta <- c(-14.26518,   42.33316,  -39.41650,  -15.53849,   46.84309,  -41.17810)
true_Sigma <- matrix(c(8.685367, -6.824736, -11.94832,
                       -6.824736, 10.382056,  14.64641,
                       -11.948320, 14.646409,  24.48155), nrow=3)

beta_length <- length(true_beta)
rand_eff_length <- 3
Nsub <- n # PASSED IN FROM THE COMMAND LINE
# For code testin, set n <- 1

seed <- sample(100:100000, 1)
set.seed(seed)

# Simulate data
rand_effects <- mvrnorm(Nsub, rep(0,3), true_Sigma)
data <- lapply(1:Nsub, function(i) {
  sub <- sample(1:Nsub, 1)
  x <- data[[sub]]$x
  z <- data[[sub]]$z
  prob <-  c(1/(1+exp(-(x %*% true_beta + z %*% rand_effects[i,]))))
  y <- rbinom(length(prob), 1, prob)
  list(x=x,z=z,y=y)
})
# Run glmer estimate
x_mat <- do.call("rbind", lapply(data, function(d) d$x))
y_mat <- unlist(lapply(data, function(d) d$y))
z_mat <- do.call("rbind", lapply(data, function(d) d$z))
sub_mat <- rep(1:Nsub, times=unlist(lapply(data, function(d) nrow(d$x))))

tt <- Sys.time()
# Note: The following command may produce a warning regarding non-convergence
# of the optimization used to estimate the GLMM.
glmer_est <- glmer(y_mat ~ x_mat[,2] + x_mat[,3] + x_mat[,4] + x_mat[,5] + x_mat[,6] + (1 + z_mat[,2] + z_mat[,3] | sub_mat), family='binomial', verbose=1)
glmer_time <- Sys.time() - tt
glmer_beta <- summary(glmer_est)$coefficients[,1]

# The way we defined the model was in terms of separate intercepts for males and females, 
# with no global intercept, so we have to  add the intercept for the coefficient for males
glmer_beta[4] <- glmer_beta[1] + glmer_beta[4]

glmer_cor <- attr(summary(glmer_est)$varcor[[1]], "correlation")
glmer_sd <- attr(summary(glmer_est)$varcor[[1]], "stddev")
glmer_Sigma <- diag(glmer_sd) %*% glmer_cor %*% diag(glmer_sd)

# Run variational EM algorithim
beta_start <- glmer_beta + rnorm(length(glmer_beta))
Sigma_start <- diag(rep(1,3))
Sigma_D_start <- log(LDL(Sigma_start)$D)
Sigma_L_start <- LDL(Sigma_start)$L
theta_start <- c(beta_start, Sigma_D_start, Sigma_L_start[lower.tri(Sigma_L_start, diag=FALSE)])


tt <- Sys.time()
vem_opt <- mixed_logit_vem(data=data, beta_start=beta_start, Sigma_start=Sigma_start, verbose=TRUE, nodes=nodes, tol=.001)
vem_time <- Sys.time() - tt

opt_par <- par_vec_to_list(c(vem_opt$beta, unlist(vem_opt$psi)), K=beta_length, P=3, n=length(data))

vem_beta <- opt_par$beta
vem_Sigma <- opt_par$Sigma
vem_Sigma_D <- log(LDL(vem_Sigma)$D)
vem_Sigma_L <- LDL(vem_Sigma)$L
vem_theta <- c(vem_beta, vem_Sigma_D, vem_Sigma_L[lower.tri(vem_Sigma_L, diag=FALSE)])

# Compute  One-step
grid <- createNIGrid(dim=3, type='GHN', level=10)
vem_score_info <- score_info(vem_beta, vem_Sigma_D, vem_Sigma_L, data=data, verbose=TRUE, nodes=nodes)

os_theta <- vem_theta + solve(vem_score_info$obs_info) %*% vem_score_info$score
os_time <- Sys.time() - tt

os_beta <- os_theta[1:6]
os_Sigma_D <- os_theta[7:9]
os_Sigma_L <- diag(rep(1,3))
os_Sigma_L[lower.tri(os_Sigma_L, diag=FALSE)] <- os_theta[10:12]
os_Sigma <- os_Sigma_L %*% diag(exp(os_Sigma_D)) %*% t(os_Sigma_L)

# Compute MLE
tt <- Sys.time()
mle_opt <- optim(theta_start, fn=log_lik_fn, gr=log_lik_grad, data=data, control=list(fnscale=-1, trace=1), method='L-BFGS-B')
mle_time <- Sys.time() - tt

mle_theta <- mle_opt$par
mle_beta <- mle_theta[1:6]
mle_Sigma_D <- mle_theta[7:9]
mle_Sigma_L <- diag(rep(1,3))
mle_Sigma_L[lower.tri(mle_Sigma_L, diag=FALSE)] <- mle_theta[10:12]
mle_Sigma <- mle_Sigma_L %*% diag(exp(mle_Sigma_D)) %*% t(mle_Sigma_L)

res <- list(seed=seed, task_id=task_id, n=n, glmer_beta=glmer_beta, vem_beta=vem_beta, mle_beta=mle_beta, os_beta=os_beta, glmer_Sigma=glmer_Sigma, vem_Sigma=vem_Sigma, mle_Sigma=mle_Sigma, os_Sigma=os_Sigma, glmer_time=glmer_time, vem_time=vem_time, os_time=os_time, mle_time=mle_time)

save(res, file=paste0('timing_results/timing_sim_n', n, '_', task_id, '.Rdata'))



