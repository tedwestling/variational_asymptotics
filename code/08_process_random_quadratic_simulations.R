###############################################################################
# REPLICATION FILES FOR:
# Beyond prediction: A framework for inference with variational approximations in mixture models
# Ted Westling and Tyler McCormick
#
# file: 08_process_random_quadratic_simulations.R
# 
# The cluster-based simulation study using 07_random_quadratic_simulations.R 
# needs to have already been run, and the 1000 results of the simulation study 
# need to have been saved in the simulation_results/ subdirectory.
#
# This script takes the output of the 1000 random quadratic simulations and produces
# a summary data frame with point estimates, confidence intervals, and the results
# of the consistency tests. It then produces Figure 2 and Table 3 of the manuscript.
#
###############################################################################

library(plyr)
library(ggplot2)
library(xtable)

source('code/mixed_logit_functions.R')

true_beta <- c(-14.26518,   42.33316,  -39.41650,  -15.53849,   46.84309,  -41.17810)
true_Sigma <- matrix(c(8.685367, -6.824736, -11.94832,
                       -6.824736, 10.382056,  14.64641,
                       -11.948320, 14.646409,  24.48155), nrow=3)
true_Sigma_D <- log(LDL(true_Sigma)$D)
true_Sigma_L <- LDL(true_Sigma)$L[lower.tri(true_Sigma, diag=FALSE)]
true_theta <- c(true_beta, true_Sigma_D, true_Sigma_L)

beta_length <- length(true_beta)
rand_eff_length <- 3
Nsub <- 8680

files <- list.files('simulation_results/',full.names = TRUE)

simulation_results <- ldply(files, function(file) {
  load(file)
  print(task_id)
  
  ##  Point estimates
  # glmer
  glmer_beta <- summary(glmer_est)$coefficients[,1]
  glmer_beta[4] <- glmer_beta[1] + glmer_beta[4]
  
  glmer_cor <- attr(summary(glmer_est)$varcor[[1]], "correlation")
  glmer_sd <- attr(summary(glmer_est)$varcor[[1]], "stddev")
  glmer_Sigma <- diag(glmer_sd) %*% glmer_cor %*% diag(glmer_sd)
  
  glmer_Sigma_D <- log(LDL(glmer_Sigma)$D)
  glmer_Sigma_L <- LDL(glmer_Sigma)$L
  glmer_theta <- c(glmer_beta, glmer_Sigma_D, glmer_Sigma_L[lower.tri(glmer_Sigma_L, diag=FALSE)])
  
  
  # vem
  opt_par <- par_vec_to_list(c(vem_opt$beta, unlist(vem_opt$psi)), K=beta_length, P=3, n=length(data))
  vem_beta <- opt_par$beta
  vem_Sigma <- opt_par$Sigma
  vem_Sigma_D <- log(LDL(vem_Sigma)$D)
  vem_Sigma_L <- LDL(vem_Sigma)$L
  vem_theta <- c(vem_beta, vem_Sigma_D, vem_Sigma_L[lower.tri(vem_Sigma_L, diag=FALSE)])
  
  # one step
  os_theta <- c(vem_theta + solve(vem_score_info$obs_info) %*% vem_score_info$score)

  ## Confidence intervals
  # glmer
  glmer_CIs <- rbind(confint(glmer_est, method='Wald')[-(1:6),], matrix(NA, nrow=6, ncol=2))
  glmer_contains <- glmer_CIs[,1] <= true_theta & true_theta <= glmer_CIs[,2] 
  
  # VEM + Sand
  vem_CIs <- cbind(vem_theta - qnorm(.975) * sqrt(diag(sand_cov)), vem_theta + qnorm(.975) * sqrt(diag(sand_cov)))
  vem_contains <- vem_CIs[,1] <= true_theta & true_theta <= vem_CIs[,2]
  
  # VEMOS + Info
  os_CIs <- try(cbind(os_theta - qnorm(.975) * sqrt(diag(solve(os_score_info$obs_info) / Nsub)), os_theta + qnorm(.975) * sqrt(diag(solve(os_score_info$obs_info) / Nsub))), silent=TRUE)
  if(class(os_CIs) == "try-error") os_CIs <- matrix(NA, nrow=12, ncol=2)
  os_contains <- os_CIs[,1] <= true_theta & true_theta <= os_CIs[,2]
 
  # Test of consistency
  marginal_pvals <- unlist(lapply(marginal_ttests, function(t) t$p.value))
  hotelling_pval <- hotelling_test$p.value
  
  data.frame(Estimate=c(glmer_theta, vem_theta, os_theta), 
             Lower_95CI=c(glmer_CIs[,1], vem_CIs[,1], os_CIs[,1]),
             Upper_95CI=c(glmer_CIs[,2], vem_CIs[,2], os_CIs[,2]),
             CI_contains=c(glmer_contains, vem_contains, os_contains),
             Parameter=rep(c("beta[0]", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "log(D[1])", "log(D[2])", "log(D[3])", "L[21]", "L[31]", "L[32]"), 3),
             Estimator=rep(c("GLMER", "VEM", "VEMOS"), each=12),
             Conf_Int=rep(c("GLMER", "VEM+Sand.", "VEMOS+Info."), each=12),
             Marginal_Pvals=c(rep(NA, 12), marginal_pvals, rep(NA, 12)),
             Hotelling_Pval=c(rep(NA, 12), rep(hotelling_pval, 12), rep(NA, 12)),
             seed=seed,
             sim=task_id)
})

save(simulation_results, file='simulation_results/random_quadratic_simulation_results.Rdata')

load('simulation_results/random_quadratic_simulation_results.Rdata')

true_beta <- c(-14.26518,   42.33316,  -39.41650,  -15.53849,   46.84309,  -41.17810)
true_Sigma <- matrix(c(8.685367, -6.824736, -11.94832,
                       -6.824736, 10.382056,  14.64641,
                       -11.948320, 14.646409,  24.48155), nrow=3)
true_Sigma_D <- log(LDL(true_Sigma)$D)
true_Sigma_L <- LDL(true_Sigma)$L[lower.tri(true_Sigma, diag=FALSE)]
true_theta <- c(true_beta, true_Sigma_D, true_Sigma_L)

true_theta_df <- data.frame(Parameter=c("beta[0]", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "log(D[1])", "log(D[2])", "log(D[3])", "L[21]", "L[31]", "L[32]"),
                            true_theta=true_theta)

levels(simulation_results$Estimator) <- c("lme4", "GVA", "GVA+OS")
(g <- ggplot(subset(simulation_results, Parameter %in% c("beta[0]", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]"))) +
  geom_boxplot(aes(Estimator, Estimate)) +
  geom_hline(data=true_theta_df[1:6,], aes(yintercept=true_theta), linetype=2) +
  facet_wrap(~Parameter, labeller=label_parsed, scales='free') +
  xlab(NULL) +
  theme_minimal())
ggsave(filename='plots/FIGURE2.png', g, width=5.5, height=3, units='in', scale=1.3)

# Compute coverages
# using na.rm = TRUE because two simulations failed to give confidence intervals
coverages <- ddply(subset(simulation_results, Parameter %in% c("beta[0]", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]")), .(Parameter), function(df){
  with(df, data.frame(MLE=mean(df$CI_contains[df$Conf_Int == "GLMER"], na.rm=TRUE),
                      GVA=mean(df$CI_contains[df$Conf_Int == "VEM+Sand."], na.rm=TRUE),
                      Onestep=mean(df$CI_contains[df$Conf_Int == "VEMOS+Info."], na.rm=TRUE)))
})

coverages2 <- data.frame(t(coverages[,2:4]))
names(coverages2) <- coverages[,1]

# MANUSCRIPT TABLE 3
xtable(coverages2)

# Summaries of p-values of marginal t-test pvalues
ddply(subset(simulation_results, Estimator=="VEM"), .(Parameter), function(df) {
  mean(df$Marginal_Pvals <.01)
})

# Summary of joint t-test p-value
summary(simulation_results$Hotelling_Pval)
