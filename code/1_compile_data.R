###############################################################################
# Variational asymptotics and sandwich estimation
#
# file: 1_compile_data.R
#
# This script takes the raw NLSY97 data and produces a data frame with the
# variables relevant for our analyses.
# The data was downloaded from https://www.nlsinfo.org/investigator/pages/search.jsp?s=NLSY97
# and saved as data/data.R. To download the exact same dataset, go to the above link 
# and upload the file data.NLSY97 (included in the replication files) as a tagset.
#
# Author: tedwestling
###############################################################################

# The data comes as a .R file that you source
setwd('data')
source('data.R')
library(MASS)

log1pex <- function(x) {
  ret <- rep(NA, length(x))
  large <- (x > 13)
  ret[large] <- x[large] 
  ret[!large] <- log1p(exp(x[!large]))
  ret
}

expit <- function(x, d=0) {
  ex <- 1/(1 + exp(-x))
  if(d==0) return(ex)
  if(d==1) return(ex - ex^2)
  if(d==2) return(ex - 3 * ex^2 + 2 * ex^3 )
  if(d==3) return(ex - 7 * ex^2 + 12 * ex^3 - 6 * ex^4)
  if(d==4) return(ex - 15 * ex^2 + 50 * ex^3 - 60 * ex^4 + 24 * ex^5)
}

logit <- function(x) log(x) - log(1-x)

library(plyr)

data_long <- ddply(new_data, .(R0000100), function(sub_data) {
  int_year <- c(1997, 1998, 1999, 2000, 2001, 2002, 2003, 2004, 2005, 2006, 2007, 2008, 2009, 2010, 2011)
  used_marijuana_sdli <- sub_data[1,c(2, 10, 13, 16, 19, 22, 26, 29, 32, 35, 38, 41, 44, 47, 50)]
  used_marijuana_last30 <- sub_data[1, c(3, 11, 14, 17, 20, 23, 27, 30, 33, 36 ,39, 42, 45, 48, 51)]
  used_marijuana_last30[used_marijuana_sdli == 0] <- 0 # Respondents not asked if used marijuana in past 30 days if they hadn't used SDLI (since last int)
  used_marijuana_last30 <- as.numeric(used_marijuana_last30 > 0)
  female_ind <- as.numeric(sub_data[1,4] == 2)
  age <- as.numeric(sub_data[1, c(7, 12, 15, 18, 21, 24, 25, 28, 31, 34, 37, 40, 43, 46, 49) ])
  age_years <- age/12
  data.frame(int_year, used_marijuana_last30, female_ind, age_years)
})

# Remove missing age, in which case subject wasn't interviewed that year
data_long <- data_long[!is.na(data_long$age_years),]

# still some reamining missing marijuana information, but not large enough to impact results in a large way
sum(is.na(data_long$used_marijuana_last30))
data_long <- data_long[!is.na(data_long$used_marijuana_last30),]

names(data_long)[1] <- 'sub'

# remove subjects with fewer than 4 time points
obs_per_sub <- table(data_long$sub)
table(obs_per_sub)
fewer <- as.numeric(names(obs_per_sub)[obs_per_sub < 4])
data_long <- subset(data_long, !(sub %in% fewer))

save(data_long, file='data_long.Rdata')


# list-ify data for random intercepts
data <- lapply(unique(data_long$sub), function(sb) {
  sub_data <- subset(data_long, sub == sb)
  y <- sub_data$used_marijuana_last30
  z <- matrix(rep(1,nrow(sub_data)), ncol=1) #cbind(1, , (sub_data$age_years/35)^2)
  age_mat <- cbind(rep(1, nrow(sub_data)), sub_data$age_years/35, (sub_data$age_years/35)^2)
  x <- cbind(sub_data$female_ind[1] * age_mat, (1-sub_data$female_ind[1]) * age_mat)
  list(x=x, z=z, y=y)
})

save(data, file='random_intercept_data.Rdata')

# list-ify data for random quadratics
data <- lapply(unique(data_long$sub), function(sb) {
  sub_data <- subset(data_long, sub == sb)
  y <- sub_data$used_marijuana_last30
  z <- cbind(1, sub_data$age_years/35, (sub_data$age_years/35)^2)
  x <- cbind(sub_data$female_ind[1] * z, (1-sub_data$female_ind[1]) * z)
  list(x=x, z=z, y=y)
})

save(data, file='random_quadratic_data.Rdata')
