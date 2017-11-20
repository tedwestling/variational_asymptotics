###############################################################################
# Variational asymptotics and sandwich estimation
#
# file: 4_process_random_intercept_simulations.R
#
# This script takes the output of the 1000 random intercept simulations and produces
# a summary data frame with point estimates, confidence intervals, and the results
# of the consistency tests. It then produces the plots and tables in the paper.
#
###############################################################################

source('code/mixed_logit_functions.R')

true_beta <- c(-10.42998,   24.72670,  -20.39093,  -12.47501,   31.81598,  -25.04588)
true_Sigma <- 4.788012
true_theta <- c(true_beta, log(true_Sigma))

beta_length <- length(true_beta)
rand_eff_length <- 1
Nsub <- 8680

files <- list.files('simulation_results/',full.names = TRUE)

simulation_results <- ldply(files, function(file) {
  load(file)
  sim <- as.numeric(regmatches(file, regexec("([0-9]+)", file))[[1]][1])
  print(sim)
  
  ##  Point estimates
  # glmer
  glmer_beta <- summary(glmer_est)$coefficients[,1]
  glmer_Sigma <- as.numeric(attr(summary(glmer_est)$varcor[[1]], "stddev"))^2
  glmer_theta <- c(glmer_beta, log(glmer_Sigma))
  
  # vem
  opt_par <- par_vec_to_list(c(vem_opt$beta, unlist(vem_opt$psi)), K=beta_length, P=1, n=length(data))
  vem_beta <- opt_par$beta
  vem_Sigma <- opt_par$Sigma
  vem_Sigma_D <- log(LDL(vem_Sigma)$D)
  vem_Sigma_L <- LDL(vem_Sigma)$L
  vem_theta <- c(vem_beta, vem_Sigma_D)
  
  # one step
  os_theta <- c(vem_theta + solve(vem_score_info$obs_info) %*% vem_score_info$score)

  ## Confidence intervals
  # glmer
  glmer_CIs <- rbind(confint(glmer_est, method='Wald')[-1,], c(NA,NA))
  glmer_contains <- glmer_CIs[,1] <= true_theta & true_theta <= glmer_CIs[,2] 
  
  # VEM + Sand
  vem_sand_CIs <- cbind(vem_theta - qnorm(.975) * sqrt(diag(sand_cov$sand_cov)), vem_theta + qnorm(.975) * sqrt(diag(sand_cov$sand_cov)))
  vem_sand_contains <- vem_sand_CIs[,1] <= true_theta & true_theta <= vem_sand_CIs[,2]
  vem_naive_CIs <- cbind(vem_theta - qnorm(.975) * sqrt(diag(-sand_cov$naive_cov)), vem_theta + qnorm(.975) * sqrt(diag(-sand_cov$naive_cov)))
  vem_naive_contains <- vem_naive_CIs[,1] <= true_theta & true_theta <= vem_naive_CIs[,2]
  
  # VEMOS + Info
  os_CIs <- cbind(os_theta - qnorm(.975) * sqrt(diag(solve(os_score_info$obs_info) / Nsub)), os_theta + qnorm(.975) * sqrt(diag(solve(os_score_info$obs_info) / Nsub)))
  os_contains <- os_CIs[,1] <= true_theta & true_theta <= os_CIs[,2]
 
  # Test of consistency
  marginal_pvals <- unlist(lapply(marginal_ttests, function(t) t$p.value))
  hotelling_pval <- hotelling_test$p.value
  hotelling_beta_pval <- hotelling_test_beta$p.value
  
  data.frame(Estimate=c(glmer_theta, vem_theta, vem_theta, os_theta), 
             Lower_95CI=c(glmer_CIs[,1], vem_sand_CIs[,1], vem_naive_CIs[,1], os_CIs[,1]),
             Upper_95CI=c(glmer_CIs[,2], vem_sand_CIs[,2], vem_naive_CIs[,2], os_CIs[,2]),
             CI_contains=c(glmer_contains, vem_sand_contains, vem_naive_contains, os_contains),
             Parameter=rep(c("beta[0]", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "log(sigma^2)"), 4),
             Estimator=rep(c("GLMER", "VEM", "VEM", "VEMOS"), each=7),
             Conf_Int=rep(c("GLMER", "VEM+Sand.", "VEM+Naive", "VEMOS+Info."), each=7),
             Marginal_Pvals=c(rep(NA,7), marginal_pvals, rep(NA, 7*2)),
             Hotelling_Pval=c(rep(NA,7), rep(hotelling_pval, 7), rep(NA, 7*2)),
             Hotelling_beta_Pval=c(rep(NA,7), rep(hotelling_beta_pval, 7), rep(NA, 7*2)),
             seed=seed,
             sim=sim)
})

joint_CI_results <- ldply(files, function(file) {
  load(file)
  sim <- as.numeric(regmatches(file, regexec("([0-9]+)", file))[[1]][1])
  print(sim)
  
  # vem
  opt_par <- par_vec_to_list(c(vem_opt$beta, unlist(vem_opt$psi)), K=beta_length, P=1, n=length(data))
  vem_beta <- opt_par$beta
  vem_Sigma <- opt_par$Sigma
  vem_Sigma_D <- log(LDL(vem_Sigma)$D)
  vem_Sigma_L <- LDL(vem_Sigma)$L
  vem_theta <- c(vem_beta, vem_Sigma_D)
  
  # one step
  os_theta <- c(vem_theta + solve(vem_score_info$obs_info) %*% vem_score_info$score)
  
  vem_sand_joint_contains <- c(as.vector(true_theta[1:2] - c(vem_theta[1:2])) %*% sand_cov$sand_cov[1:2,1:2] %*% (true_theta[1:2] - c(vem_theta[1:2])) <= qchisq(.95, 2))
  vem_naive_joint_contains <- c(as.vector(true_theta[1:2] - c(vem_theta[1:2])) %*% (-sand_cov$naive_cov)[1:2,1:2] %*% (true_theta[1:2] - c(vem_theta[1:2])) <= qchisq(.95, 2))
  os_info_joint_contains <- c(as.vector(true_theta[1:2] - c(os_theta[1:2])) %*% (solve(os_score_info$obs_info) / Nsub)[1:2,1:2] %*% (true_theta[1:2] - c(os_theta[1:2])) <= qchisq(.95, 2))
  data.frame(contains=c(vem_sand_joint_contains, vem_naive_joint_contains, os_info_joint_contains), method=c("vem_sand", "vem_naive", "os_info"))
})
  

save(simulation_results, file='random_intercept_simulation_results.Rdata')

## 
source('code/mixed_logit_functions.R')

true_beta <- c(-10.42998,   24.72670,  -20.39093,  -12.47501,   31.81598,  -25.04588)
true_Sigma <- 4.788012
true_theta <- c(true_beta, log(true_Sigma))

beta_length <- length(true_beta)
rand_eff_length <- 1
Nsub <- 8680

load('simulation_results/random_intercept_simulation_results.Rdata')

true_theta_df <- data.frame(Parameter=c("beta[0]", "beta[1]", "beta[2]", "beta[3]", "beta[4]", "beta[5]", "log(sigma^2)"),
                            true_theta=true_theta)

library(ggplot2)
(g <- ggplot(simulation_results) +
  geom_boxplot(aes(Estimator, Estimate)) +
  geom_hline(data=true_theta_df, aes(yintercept=true_theta), linetype=2) +
  facet_wrap(~Parameter, labeller=label_parsed, scales='free') +
  xlab(NULL) +
  scale_x_discrete(breaks=c("GLMER", "VEM", "VEMOS"), labels=c("lme4", "GVA", "GVA + OS")) +
  theme_minimal())
ggsave(filename='plots/random_intercept_sim_estimates.png', g, width=5.5, height=3, units='in', scale=1.3)

# Compute MSEs
mses <- ddply(simulation_results, .(Parameter), function(df) {
  truth <- true_theta_df$true_theta[true_theta_df$Parameter == df$Parameter[1]]
  with(df, data.frame(GLMER_MSE=mean((Estimate[Estimator == "GLMER"] - truth)^2),
                      VEM_MSE=mean((Estimate[Estimator == "VEM"] - truth)^2),
                      OS_MSE=mean((Estimate[Estimator == "VEMOS"] - truth)^2)))
})

variances <- ddply(simulation_results, .(Parameter), function(df) {
  truth <- true_theta_df$true_theta[true_theta_df$Parameter == df$Parameter[1]]
  with(df, data.frame(GLMER_VAR=var(Estimate[Estimator == "GLMER"]),
                      VEM_VAR=var(Estimate[Estimator == "VEM"]),
                      OS_VAR=var(Estimate[Estimator == "VEMOS"])))
})

variances2 <- data.frame(t(variances[,-1]))
names(variances2) <- variances[,1]
xtable(variances2)


coverages <- ddply(subset(simulation_results, Parameter != "log(sigma^2)" | Estimator != "GLMER"), .(Parameter, Conf_Int), function(df){
  ptest <- prop.test(x=sum(df$CI_contains), n=nrow(df))
  data.frame(Coverage=ptest$estimate, Lower_95CI=ptest$conf.int[1], Upper_95CI=ptest$conf.int[2])
})


(g <- ggplot(subset(coverages, Parameter != "log(sigma^2)")) + 
  geom_point(aes(Conf_Int, Coverage)) +
  geom_errorbar(aes(Conf_Int, ymin=Lower_95CI, ymax=Upper_95CI)) +
  geom_hline(yintercept=.95, linetype=2) +
  facet_wrap(~Parameter, nrow=2, labeller=label_parsed) +
  ylim(c(.85, 1)) +
  xlab(NULL) +
  theme_minimal())
ggsave(filename='plots/random_intercept_sim_CIs.png', g, width=5.5, height=3, units='in', scale=1.3)

# histograms of p-values of marginal t-test pvalues
ggplot(subset(simulation_results, Estimator=="GLMER")) + 
  geom_histogram(aes(Marginal_Pvals, y=..density..)) + 
  facet_wrap(~Parameter, labeller=label_parsed, scales='free') +
  theme_minimal()
ddply(subset(simulation_results, Estimator=="GLMER"), .(Parameter), function(df) {
  mean(df$Marginal_Pvals <.01)
})

ggplot(subset(simulation_results, Estimator=="GLMER" & Parameter == "beta[0]")) + 
  geom_histogram(aes(Hotelling_Pval, y=..density..)) + 
  theme_minimal()
summary(simulation_results$Hotelling_Pval)
