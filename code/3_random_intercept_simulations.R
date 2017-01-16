###############################################################################
# Variational asymptotics and sandwich estimation
#
# file: 3_random_intercept_simulations.R
#
# This script runs a single logistic regression with random intercepts simulation.
# The true parameters were drawn from the estimates for the NLSY97 data.
# MLE (with glmer), variational EM, and variational EM + one step are used
# to estimate the parameters, and the asymptotic covariance matrices are calculated
# for ecah. Consistency of the variational EM algorithm is also assessed.
# In addition to running the simulations in parallel on a cluster, 12 cores are used
# per simulation to speed up computation. Here is the slurm code used to run all 1000 
# simulations:
# !/bin/bash
# #SBATCH --job-name random_intercept_sim      # Set a name for your job. This is especially useful if you have multiple jobs queued.
# #SBATCH --partition medium     # Slurm partition to use
# #SBATCH --ntasks 1         # Number of tasks to run. By default, one CPU core will be allocated per task
# #SBATCH --cpus-per-task=12
# #SBATCH --time 2-12:00        # Wall time limit in D-HH:MM
# #SBATCH --mem-per-cpu=800     # Memory limit for each tasks (in MB)
# #SBATCH -o outs/rand_int_sim_%A_%a.out    # File to which STDOUT will be written
# #SBATCH -e outs/rand_int_sim_%A_%a.err    # File to which STDERR will be written
# #SBATCH --mail-type=ALL       # Type of email notification- BEGIN,END,FAIL,ALL
# #SBATCH --mail-user=tgwest@uw.edu # Email to which notifications will be sent
# 
# R CMD BATCH --vanilla --no-restore --no-save /code/3_random_intercept_simulations.R  outs/random_intercept_sim_"${SLURM_ARRAY_TASK_ID}".Rout
#
# Author: tedwestling
###############################################################################

rm(list=ls())
setwd('/homes/tgwest/random_intercepts')
source('code/mixed_logit_functions.R')
nodes <- 12 # Run on the cluster
task_id <- as.numeric(Sys.getenv('SLURM_ARRAY_TASK_ID'))

# Load compiled data
load('data/random_intercept_data.Rdata')

true_beta <- c(-10.42998,   24.72670,  -20.39093,  -12.47501,   31.81598,  -25.04588)
true_Sigma <- 4.788012

beta_length <- length(true_beta)
rand_eff_length <- 1
Nsub <- length(data)

seed <- sample(100:100000, 1)
print(seed)
set.seed(seed)

# Simulate data
rand_intercepts <- rnorm(Nsub, 0, sqrt(true_Sigma))
data <- lapply(1:Nsub, function(i) {
  sub <- sample(1:Nsub, 1)
  x <- data[[sub]]$x
  z <- data[[sub]]$z
  prob <-  c(1/(1+exp(-(x %*% true_beta + z * rand_intercepts[i]))))
  y <- rbinom(length(prob), 1, prob)
  list(x=x,z=z,y=y)
})

# Run glmer estimate
x_mat <- do.call("rbind", lapply(data, function(d) d$x))
y_mat <- unlist(lapply(data, function(d) d$y))
z_mat <- do.call("rbind", lapply(data, function(d) d$z))
sub_mat <- rep(1:Nsub, times=unlist(lapply(data, function(d) nrow(d$x))))

glmer_est <- glmer(y_mat ~ 0 + x_mat[,1] + x_mat[,2] + x_mat[,3] + x_mat[,4] + x_mat[,5] + x_mat[,6] + (1 | sub_mat), family='binomial', verbose=1)
glmer_beta <- summary(glmer_est)$coefficients[,1]

glmer_Sigma <- as.numeric(attr(summary(glmer_est)$varcor[[1]], "stddev"))^2
glmer_theta <- c(glmer_beta, log(glmer_Sigma))

# Run variational EM algorithim
beta_start <- glmer_beta + rnorm(length(glmer_beta))
Sigma_start <- exp(log(glmer_Sigma) + rnorm(1, 0, .4))

psi_hats <- mclapply(1:Nsub, function(i) {
  if(i %% 500 == 0) cat("Obs", i, "of", Nsub, "\n")
  one_data <- data[[i]]
  
  this_opt <- optim(rep(0, 2), fn=one_elbo_theta0, gr=one_elbo_grad_theta0, beta=beta_start, Sigma=matrix(Sigma_start, nrow=1,ncol=1), Yi=one_data$y, Zi=one_data$z, Xi=one_data$x, control=list(fnscale=-1), method='BFGS')
  mi <- this_opt$par[1]
  Di <- this_opt$par[2]
  Li <- diag(1)
  
  list(mi=mi, Di=Di, Li=Li)
}, mc.cores=nodes)

start_m <- matrix(sapply(1:Nsub, function(i) psi_hats[[i]]$mi),ncol=1)
start_D <- matrix(sapply(1:Nsub, function(i) psi_hats[[i]]$Di),ncol=1)
start_L <- matrix(lapply(1:Nsub, function(i) psi_hats[[i]]$Li),ncol=1)

startvec <- par_list_to_vec(beta=beta_start, m=start_m, D=start_D, L=start_L)
start_par <- par_vec_to_list(startvec, K=6, P=1, n=Nsub)

vem_opt <- mixed_logit_vem(data=data, beta_start=beta_start, Sigma_start=start_par$Sigma, verbose=FALSE, nodes=nodes)

opt_par <- par_vec_to_list(c(vem_opt$beta, unlist(vem_opt$psi)), K=beta_length, P=1, n=length(data))

vem_beta <- opt_par$beta
vem_Sigma <- opt_par$Sigma
vem_Sigma_D <- log(LDL(vem_Sigma)$D)
vem_Sigma_L <- LDL(vem_Sigma)$L
vem_theta <- c(vem_beta, vem_Sigma_D, vem_Sigma_L[lower.tri(vem_Sigma_L, diag=FALSE)])

# Compute sandwich covariance
sand_cov <- sandwich_cov(beta=vem_beta, Sigma_D=vem_Sigma_D, Sigma_L=vem_Sigma_L, m=opt_par$m, D=opt_par$D, L=opt_par$L, data=data, verbose=TRUE, nodes=nodes)

# One-step correction
grid <- createNIGrid(dim=1, type='GHN', level=10)
vem_score_info <- score_info(vem_beta, vem_Sigma_D, vem_Sigma_L, data=data, verbose=TRUE, nodes=nodes)

os_theta <- vem_theta + solve(vem_score_info$obs_info) %*% vem_score_info$score
os_beta <- os_theta[1:length(glmer_beta)]
os_Sigma_D <- os_theta[length(glmer_beta)+1]
os_Sigma <- exp(os_Sigma_D)
os_Sigma_L <- 1

os_score_info <- score_info(os_beta, os_Sigma_D, os_Sigma_L, data=data, verbose=TRUE, nodes=nodes)

n_sim <- 100000
sub_effs <- rnorm(n_sim, 0, sqrt(vem_Sigma))

library(parallel)
start_time <- Sys.time()
prof_func_hats <- mclapply(1:n_sim, function(i) {
  if(i %% 500 == 0) {
    time_left <- (c(difftime(Sys.time(), start_time, units="secs")) / i) * (n_sim - i)
    cat("Obs", i, "of", n_sim, ",", round(c(time_left)/60, 2), "minutes left\n")
  }
  samp <- sample(1:length(data), 1)
  z <- data[[samp]]$z
  x <- data[[samp]]$x
  prob <-  c(1/(1+exp(-(x %*% vem_beta + z * sub_effs[i]))))
  y <- rbinom(length(prob), 1, prob)
  this_opt <- optim(rep(0, 2), fn=one_elbo_theta0, gr=one_elbo_grad_theta0, beta=vem_beta, Sigma=vem_Sigma, Yi=y, Zi=z, Xi=x, control=list(fnscale=-1), method='L-BFGS-B')
  
  mi <- this_opt$par[1]
  Di <- this_opt$par[2]
  Li <- matrix(1, nrow=1)

  grad <- one_elbo_thetagrad(beta=vem_beta, Sigma=vem_Sigma, mi=mi, Di=Di, Li=Li, Yi=y, Zi=z, Xi=x)
  list(this_opt, grad)
}, mc.cores=nodes)

prof_func_hats_mat <- matrix(unlist(lapply(prof_func_hats, function(l) l[[2]])), ncol=beta_length+1, byrow=TRUE)
marginal_ttests <- apply(prof_func_hats_mat, 2, t.test)
library(DescTools)
hotelling_test <- HotellingsT2Test(prof_func_hats_mat)

save(seed, glmer_est, vem_opt, sand_cov, vem_score_info, os_score_info, hotelling_test, marginal_ttests, file=paste0('simulation_results/rand_intercepts_sim', task_id, '.Rdata'))