###############################################################################
# REPLICATION FILES FOR:
# Beyond prediction: A framework for inference with variational approximations in mixture models
# Ted Westling and Tyler McCormick
#
# file: 04_make_NLSY_plots.R
# 
# The files 02_run_random_intercepts.R and 03_run_random_quadratics.R
# need to have been run, which will have created the files 
# estimates/rand_intercepts_estimate.Rdata and 
# estimates/rand_quadratics_estimate.Rdata
#
# This script produces Figure 4 in the manuscript and additional results of the
# data analysis described in the text.
#
# Author: Ted Westling (ted.westling@gmail.com)
###############################################################################

library(ggplot2)
library(DescTools)

source('code/mixed_logit_functions.R')
beta_length <- 6
Nsub <- 8680

# The file 02_run_random_intercepts needs to have been run, which will have created the file 
# estimates/rand_intercepts_estimate.Rdata
load('estimates/rand_intercepts_estimate.Rdata')

opt_par <- par_vec_to_list(c(vem_opt$beta, unlist(vem_opt$psi)), K=beta_length, P=1, n=Nsub)

vem_beta <- opt_par$beta
vem_Sigma <- opt_par$Sigma
vem_Sigma_D <- log(LDL(vem_Sigma)$D)
vem_Sigma_L <- LDL(vem_Sigma)$L
vem_theta <- c(vem_beta, vem_Sigma_D, vem_Sigma_L[lower.tri(vem_Sigma_L, diag=FALSE)])


int_os_theta <- vem_theta + solve(vem_score_info$obs_info) %*% vem_score_info$score
int_os_beta <- int_os_theta[1:6]
int_os_cov <- solve(os_score_info$obs_info) / Nsub
int_os_CIs <- cbind(int_os_theta - qnorm(.975) * sqrt(diag(int_os_cov)), int_os_theta + qnorm(.975) * sqrt(diag(int_os_cov)))
int_os_CIs <- int_os_CIs[1:6,]

# The file 04_run_random_quadratics.R needs to have been run, which will have created the file 
# estimates/rand_quadratics_estimate.Rdata
load('estimates/rand_quadratics_estimate.Rdata')
opt_par <- par_vec_to_list(c(vem_opt$beta, unlist(vem_opt$psi)), K=beta_length, P=3, n=Nsub)

vem_beta <- opt_par$beta
vem_Sigma <- opt_par$Sigma
vem_Sigma_D <- log(LDL(vem_Sigma)$D)
vem_Sigma_L <- LDL(vem_Sigma)$L
vem_theta <- c(vem_beta, vem_Sigma_D, vem_Sigma_L[lower.tri(vem_Sigma_L, diag=FALSE)])

quad_os_theta <- c(vem_theta + solve(vem_score_info$obs_info) %*% vem_score_info$score)
quad_os_beta <- quad_os_theta[1:6]
quad_os_cov <- solve(os_score_info$obs_info) / Nsub
quad_os_CIs <- cbind(quad_os_theta - qnorm(.975) * sqrt(diag(quad_os_cov)), quad_os_theta + qnorm(.975) * sqrt(diag(quad_os_cov)))
quad_os_CIs <- quad_os_CIs[1:6,]

full_ests <- data.frame(estimate=c(int_os_beta, quad_os_beta),
                        lower=c(int_os_CIs[,1], quad_os_CIs[,1]),
                        upper=c(int_os_CIs[,2], quad_os_CIs[,2]),
                        model=rep(c("Rand. int.", "Rand. quad."), each=6),
                        parameter=rep(c("beta[0]", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]"), 2))

age_grid <- seq(10, 35, length.out = 100)
predicted_curve <- function(ages, beta) 1/(1+exp(-(beta[1] + beta[2] * (ages/35) + beta[3] * (ages/35)^2)))
pointwise_band <- function(x, est, cov, p) {
  sds <- sapply(x, function(xx) {
    est_curve <- est[1] + est[2] * (xx/35) + est[3] * (xx/35)^2
    sd <- sqrt(t(c(1, xx/35, (xx/35)^2)) %*% cov %*% c(1, xx/35, (xx/35)^2))
    expit(qnorm(p, est_curve, sd))
  })
}

predicted <- data.frame(Age=rep(age_grid, 2),
                        pred=c(predicted_curve(age_grid, int_os_beta[1:3]), 
                               predicted_curve(age_grid, int_os_beta[4:6]),
                               predicted_curve(age_grid, quad_os_beta[1:3]), 
                               predicted_curve(age_grid, quad_os_beta[4:6])),
                        lower=c(pointwise_band(age_grid, int_os_beta[1:3], int_os_cov[1:3, 1:3], .025),  
                                pointwise_band(age_grid, int_os_beta[4:6], int_os_cov[4:6, 4:6], .025),
                                pointwise_band(age_grid, quad_os_beta[1:3], quad_os_cov[1:3, 1:3], .025),  
                                pointwise_band(age_grid, quad_os_beta[4:6], quad_os_cov[4:6, 4:6], .025)),
                        upper=c(pointwise_band(age_grid, int_os_beta[1:3], int_os_cov[1:3, 1:3], .975),  
                                pointwise_band(age_grid, int_os_beta[4:6], int_os_cov[4:6, 4:6], .975),
                                pointwise_band(age_grid, quad_os_beta[1:3], quad_os_cov[1:3, 1:3], .975),  
                                pointwise_band(age_grid, quad_os_beta[4:6], quad_os_cov[4:6, 4:6], .975)),
                        Sex=rep(c("Female", "Male", "Female", "Male"), each=length(age_grid)),
                        Model=rep(c("Rand. int.", "Rand. quad."), each=2*length(age_grid)))

# Make Figure 4 (predicted probability of marijuana use by age for the average Male/Female from the two models)
(g <- ggplot(predicted) +
    geom_ribbon(aes(Age, ymin=lower, ymax=upper, group=Model), alpha=.5) +
    geom_line(aes(Age, pred,  linetype=Model)) +
    facet_grid(~Sex) + 
    theme_minimal() +
    ylab("Predicted probability"))
ggsave(filename='plots/FIGURE4.png', g, width=6.5, height=3, units='in', scale=1.3)

# Estimate peak usage times and 95% intervals using the delta method
peak <- function(beta, cov, level=.95) {
  est <- -35 * beta[2] / (2 * beta[3])
  sd <- sqrt(t(c(0, -35/(2 * beta[3]),-35 * beta[2] / (2 * beta[3]^2))) %*% cov %*% c(0, -35/(2 * beta[3]),-35 * beta[2] / (2 * beta[3]^2)))
  list(estimate=est, conf.int=c(est - qnorm(1-(1-level)/2) * sd, est + qnorm(1-(1-level)/2) * sd))
}
peak(int_os_beta[1:3], int_os_cov[1:3, 1:3])
peak(int_os_beta[4:6], int_os_cov[4:6, 4:6])

peak(quad_os_beta[1:3], quad_os_cov[1:3, 1:3])
peak(quad_os_beta[4:6], quad_os_cov[4:6, 4:6])

# Empirical tests of consistency

prof_func_hats_mat <- matrix(unlist(lapply(prof_func_hats, function(l) l[[2]])), ncol=12, byrow=TRUE)
apply(prof_func_hats_mat, 2, t.test)

HotellingsT2Test(prof_func_hats_mat)

# An additional figure (not included in the article) illustrating the point estimates and CIs for each fixed effect parameter from the two models

(g <- ggplot(full_ests) +
    geom_point(aes(model, estimate)) + 
    geom_errorbar(aes(model, ymin=lower, ymax=upper)) + 
    theme_minimal() +
    facet_wrap(~parameter, scales='free', labeller=label_parsed) +
    ylab("Estimate + 95% CI") +
    xlab(NULL))

