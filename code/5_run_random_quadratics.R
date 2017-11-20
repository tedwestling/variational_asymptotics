###############################################################################
# Variational asymptotics and sandwich estimation
#
# file: 5_run_random_quadratics.R
#
# This script computes the various estimates of the logistic random quadratic model
# for the marijuana use variables of the NLSY97 data. There are six fixed
# effects corresponding to the intercept, slope, and quadratic effect of time
# for the average male and female separately.
#
###############################################################################

source('code/mixed_logit_functions.R')
nodes <- 1 # we ran this file on a personal computer with no parallelization

# Load compiled data
load('data/random_quadratic_data.Rdata')

beta_length <- ncol(data[[1]]$x)
rand_eff_length <- 3
Nsub <- length(data)

# Run glmer estimate
x_mat <- do.call("rbind", lapply(data, function(d) d$x))
y_mat <- unlist(lapply(data, function(d) d$y))
z_mat <- do.call("rbind", lapply(data, function(d) d$z))
sub_mat <- rep(1:Nsub, times=unlist(lapply(data, function(d) nrow(d$x))))

glmer_est <- glmer(y_mat ~ 0 + x_mat[,1] + x_mat[,2] + x_mat[,3] + x_mat[,4] + x_mat[,5] + x_mat[,6] + (1 + z_mat[,2] + z_mat[,3] | sub_mat), family='binomial', verbose=1)
glmer_beta <- summary(glmer_est)$coefficients[,1]

glmer_cor <- attr(summary(glmer_est)$varcor[[1]], "correlation")
glmer_sd <- attr(summary(glmer_est)$varcor[[1]], "stddev")
glmer_Sigma <- diag(glmer_sd) %*% glmer_cor %*% diag(glmer_sd)

glmer_Sigma_D <- log(LDL(glmer_Sigma)$D)
glmer_Sigma_L <- LDL(glmer_Sigma)$L
glmer_theta <- c(glmer_beta, glmer_Sigma_D, glmer_Sigma_L[lower.tri(glmer_Sigma_L, diag=FALSE)])

# Run variational EM algorithim
set.seed(187112)
beta_start <- glmer_beta + rnorm(length(glmer_beta))
Sigma_start <- diag(rep(1,3))

vem_opt <- mixed_logit_vem(data=data, beta_start=glmer_beta, Sigma_start=diag(rep(1,3)), verbose=TRUE, nodes=nodes)

opt_par <- par_vec_to_list(c(vem_opt$beta, unlist(vem_opt$psi)), K=beta_length, P=3, n=length(data))

vem_beta <- opt_par$beta
vem_Sigma <- opt_par$Sigma
vem_Sigma_D <- log(LDL(vem_Sigma)$D)
vem_Sigma_L <- LDL(vem_Sigma)$L
vem_theta <- c(vem_beta, vem_Sigma_D, vem_Sigma_L[lower.tri(vem_Sigma_L, diag=FALSE)])

# Compute sandwich covariance
sand_cov <- sandwich_cov(beta=vem_beta, Sigma_D=vem_Sigma_D, Sigma_L=vem_Sigma_L, m=opt_par$m, D=opt_par$D, L=opt_par$L, data=data, verbose=TRUE, nodes=nodes)

# Compute 
grid <- createNIGrid(dim=3, type='GHN', level=10)
vem_score_info <- score_info(vem_beta, vem_Sigma_D, vem_Sigma_L, data=data, verbose=TRUE, nodes=nodes)

os_theta <- vem_theta + solve(vem_score_info$obs_info) %*% vem_score_info$score
os_beta <- os_theta[1:6]
os_Sigma_D <- os_theta[7:9]
os_Sigma_L <- diag(rep(1,3))
os_Sigma_L[lower.tri(os_Sigma_L, diag=FALSE)] <- os_theta[10:12]
os_Sigma <- os_Sigma_L %*% diag(exp(os_Sigma_D)) %*% t(os_Sigma_L)

os_score_info <- score_info(os_beta, os_Sigma_D, os_Sigma_L, data=data, verbose=TRUE, nodes=nodes)

n_sim <- 100000
set.seed(99099092)
sub_effs <- mvrnorm(n_sim, rep(0,3), vem_Sigma)

start_time <- Sys.time()
prof_func_hats <- mclapply(1:n_sim, function(i) {
  if(i %% 500 == 0) {
    #time_left <- (c(difftime(Sys.time(), start_time, units="secs")) / i) * (n_sim - i)
    cat("Obs", i, "of", n_sim)#, ",", round(c(time_left)/60, 2), "minutes left\n")
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

prof_func_hats_mat <- matrix(unlist(lapply(prof_func_hats, function(l) l[[2]])), ncol=12, byrow=TRUE)
apply(prof_func_hats_mat, 2, t.test)
library(DescTools)
HotellingsT2Test(prof_func_hats_mat)

save(glmer_est, vem_opt, sand_cov, vem_score_info, os_score_info, prof_func_hats, file='estimates/rand_quadratics_estimate.Rdata')

######## MAKE PLOTS ########

load('estimates/rand_intercepts_estimate.Rdata')
opt_par <- par_vec_to_list(c(vem_opt$beta, unlist(vem_opt$psi)), K=beta_length, P=1, n=length(data))

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



load('estimates/rand_quadratics_estimate.Rdata')
opt_par <- par_vec_to_list(c(vem_opt$beta, unlist(vem_opt$psi)), K=beta_length, P=3, n=length(data))

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

library(ggplot2)
(g <- ggplot(full_ests) +
  geom_point(aes(model, estimate)) + 
  geom_errorbar(aes(model, ymin=lower, ymax=upper)) + 
  theme_minimal() +
  facet_wrap(~parameter, scales='free', labeller=label_parsed) +
  ylab("Estimate + 95% CI") +
  xlab(NULL))
ggsave(filename='plots/both_nlsy97_ests.png', g, width=6.5, height=4, units='in', scale=1.3)


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

(g <- ggplot(predicted) +
   geom_ribbon(aes(Age, ymin=lower, ymax=upper, group=Model), alpha=.5) +
   geom_line(aes(Age, pred,  linetype=Model)) +
   facet_grid(~Sex) + 
   theme_minimal() +
   ylab("Predicted probability"))
ggsave(filename='plots/both_nlsy97_pred.png', g, width=6.5, height=3, units='in', scale=1.3)

peak <- function(beta, cov, level=.95) {
  est <- -35 * beta[2] / (2 * beta[3])
  sd <- sqrt(t(c(0, -35/(2 * beta[3]),-35 * beta[2] / (2 * beta[3]^2))) %*% cov %*% c(0, -35/(2 * beta[3]),-35 * beta[2] / (2 * beta[3]^2)))
  list(estimate=est, conf.int=c(est - qnorm(1-(1-level)/2) * sd, est + qnorm(1-(1-level)/2) * sd))
}
peak(int_os_beta[1:3], int_os_cov[1:3, 1:3])
peak(int_os_beta[4:6], int_os_cov[4:6, 4:6])

peak(quad_os_beta[1:3], quad_os_cov[1:3, 1:3])
peak(quad_os_beta[4:6], quad_os_cov[4:6, 4:6])
