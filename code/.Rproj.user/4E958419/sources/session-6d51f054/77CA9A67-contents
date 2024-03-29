# ==============================================================================
# FInancaial Economitrics Coursework
# ==============================================================================
#
# Author: 10710007
# Version: 13-03-2024
#
# ==============================================================================

# Packages
library(lubridate)
# Data
rm(list = ls())
setwd("/home/oddish3/Documents/R_folder/MSc/FE/FE-coursework")
data <- read.csv("./data/group_11.csv")
source("../FE-utils.R")

# Script
# ==============================================================================
data$date <- as.Date(data$date)
plot(data$log.return, type = "l", col = "darkgreen")
abline(h = 0, v = 250, col = "red")
# dev.off()


# ------------------------------------------
#           Practical Exercise
# ------------------------------------------

# Use the first 250 observations as an in-sample (estimation) period and the last 250 observations as out of sample forecasting
in_sample <- as.matrix(data[1:250, 2])
out_sample <- as.matrix(data[(251:nrow(data)), 2])

# a) investigating statistical properties of the in-sample data
# i) descriptive stats
a1_results <- dstats(in_sample)

# ii) moment tests
a2_results <- test_moment(in_sample)

# lbq test
in_sample2 <- in_sample^2
lbq1 <- Box.test(in_sample, lag = 21, type = "Ljung-Box", fitdf = p + q + 2)
lbq1$p.value

lbq2 <- Box.test(in_sample^2, lag = 21, type = "Ljung-Box", fitdf = p + q + 2) # fitdf is the number of parameters estimated ???
lbq2$p.value

# Plot SACF and SPACF
par(mfrow = c(1, 2))
Acf(in_sample, lag.max = 25)
Pacf(in_sample, lag.max = 25)

# Assume that the conditional mean of the return series is constant: rt = c + εt , εt = σt zt . Define
#the residual series as et = rt − ĉ. Use the et series to estimate the following conditional variance

# b) estimating the conditional variance
# i) GARCH(1,1)
library(tseries)
b1_results <- garch(in_sample, order = c(1, 1))

