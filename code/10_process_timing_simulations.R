###############################################################################
# REPLICATION FILES FOR:
# Beyond prediction: A framework for inference with variational approximations in mixture models
# Ted Westling and Tyler McCormick
#
# file: 10_process_timing_simulations.R
#
# The 300 timing results should have already been run using 09_timing_simulations.R, 
# and the results saved in the directory timing_results/
#
# This file processes the timing results and makes Figure 3 found in the 
# manuscript.
#
# Author: Ted Westling
###############################################################################

library(plyr)
library(ggplot2)
library(reshape)

# The timing results should have already been run using 8_timing_simulations.R, 
# and the results saved in the directory timing_results/.
files <- list.files('timing_results',full.names = TRUE)

results <- ldply(files, function(file) {
  load(file)
  data.frame(n=res$n,
             id=res$task,
             glmer_time=as.numeric(res$glmer_time, units='mins'),
             vem_time=as.numeric(res$vem_time, units='mins'),
             os_time=as.numeric(res$os_time, units='mins'),
             mle_time=as.numeric(res$mle_time, units='mins'))
})
save(results, file='simulation_results/timing_results.Rdata')

load('simulation_results/timing_results.Rdata')

res2 <- melt(results, id.vars='n', measure.vars = c("glmer_time", "vem_time", "os_time", "mle_time"))
levels(res2$variable) <- c("lme4", "GVA", "GVA+OS", "MLE")

(g <- ggplot(res2) +
    geom_boxplot(aes(as.factor(n), value)) +
    facet_grid(~variable) +
    ylab("Time (minutes)") +
    theme_minimal() +
    coord_cartesian(ylim=c(0,80)) +
    xlab("n"))
ggsave("plots/FIGURE3.png", width=5, height=2, units='in', scale=1.3)