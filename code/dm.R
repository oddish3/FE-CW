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
library(forecast)
library(rugarch)
# Data
rm(list = ls())
# setwd("/home/oddish3/Documents/R_folder/MSc/FE/FE-coursework")
data <- read.csv("../data/group_11.csv")
source("fineco_fun.R")

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

# a) investigating statistical properties of the in-sample data ----
# i) descriptive stats
a1_results <- dstats(in_sample)

# ii-v) moment tests
a2_results <- test_moment(in_sample)

# vi) lbq test
lbq1 <- Box.test(in_sample, lag = 21, type = "Ljung-Box", fitdf = 0)
lbq1$p.value

lbq2 <- Box.test(in_sample^2, lag = 21, type = "Ljung-Box", fitdf = 0) # fitdf is the number of parameters estimated ???
lbq2$p.value

# vii) Plot SACF and SPACF
par(mfrow = c(1, 2))
Acf(in_sample, lag.max = 25)
Pacf(in_sample, lag.max = 25)

# Assume that the conditional mean of the return series is constant: rt = c + εt , εt = σt zt . Define
#the residual series as et = rt − ĉ. Use the et series to estimate the following conditional variance

# b) estimating the conditional variance -----
# i) GARCH(1,1) with zt ~ N(0,1)
spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                   distribution.model = "norm")
fit1 <- ugarchfit(spec = spec, data = in_sample)

estimates1 <- fit1@fit$robust.matcoef[,1]  # This extracts the "Estimate" column
p_values1<- fit1@fit$robust.matcoef[,4]  # This extracts the "Pr(>|t|)" column


# GARCH(1,1) with zt ~ tv
spec <- ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                   distribution.model = "std")  # 'std' for Student's t-distribution
fit2 <- ugarchfit(spec = spec, data = in_sample)

estimates2 <- fit2@fit$robust.matcoef[,1]  # This extracts the "Estimate" column
p_values2 <- fit2@fit$robust.matcoef[,4]  # This extracts the "Pr(>|t|)" column

# ii) GJR-GARCH (1, 1) with zt ~ N(0,1) 
spec<- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = TRUE),
  distribution.model = "norm" # Standard normal distribution for innovations
)

fit3 <- ugarchfit(spec = spec, data = in_sample)
estimates3 <- fit3@fit$robust.matcoef[,1]  # This extracts the "Estimate" column
p_values3 <- fit3@fit$robust.matcoef[,4]  # This extracts the "Pr(>|t|)" column

# GJR-MODEL (1, 1) with zt ~ tv
spec <- ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
  distribution.model = "std"  # 'std' for Student's t-distribution
)
fit4 <- ugarchfit(spec = spec, data = in_sample)
estimates4 <- fit4@fit$robust.matcoef[,1]  # This extracts the "Estimate" column
p_values4 <- fit4@fit$robust.matcoef[,4] 

# c) plotting NIC ----
# For GARCH(1,1) with zt ~ N(0,1)
w1 <- estimates1["omega"]
a1 <- estimates1["alpha1"]
b1 <- estimates1["beta1"]

# For GARCH(1,1) with zt ~ tv
w2 <- estimates2["omega"]
a2 <- estimates2["alpha1"]
b1 <- estimates2["beta1"]
v <- estimates2["shape"]

# For GJR-GARCH(1,1) with zt ~ N(0,1)
w3 <- estimates3["omega"]
a3 <- estimates3["alpha1"]
b3 <- estimates3["beta1"]
g3 <- estimates3["gamma1"]

# For GJR-GARCH(1,1) with zt ~ tv
w4 <- estimates4["omega"]
a4 <- estimates4["alpha1"]
b4 <- estimates4["beta1"]
g4 <- estimates4["gamma1"]
v4 <- estimates4["shape"]

# Unconditional variance from the first GARCH(1,1) model
ve1 <- w1 / (1 - a1 - b1)

# NIC
T = 500
e = seq(-5, 5, length.out = T)  # Grid of shocks epsilon
nicG1 = nicG2 = nicGJR1 = nicGJR2 = rep(0, T)  # Initialize NIC for each model

# Calculate NIC for each model
for (t in 1:T) {
  nicG1[t] = w1 + b1 * ve1 + a1 * e[t]^2  # GARCH(1,1) with zt ~ N(0,1)
  nicG2[t] = w2 + b1 * ve1 + a2 * e[t]^2  # GARCH(1,1) with zt ~ tv

  if (e[t]>0){
    nicGJR1[t] = w3 + b3 * ve1 + a3 * e[t]^2 # GJR-GARCH(1,1) with zt ~ N(0,1)
    nicGJR2[t] = w4 + b4 * ve1 + a4 * e[t]^2 # GJR-GARCH(1,1) with zt ~ tv
  }
  else{
    nicGJR1[t] = w3 + b3 * ve1 + (a3+g3) * e[t]^2 # GJR-GARCH(1,1) with zt ~ N(0,1)
    nicGJR2[t] = w4 + b4 * ve1 + (a4+g4) * e[t]^2 # GJR-GARCH(1,1) with zt ~ tv
  }
}

# d) best model by analysing standarised residuals ----

# Extract standardized residuals from each model
std_resid_gjr_norm <- fit1@fit[["residuals"]] / sqrt(fit1@fit[["sigma"]] )
std_resid_gjr_t <- fit2@fit[["residuals"]] / sqrt(fit2@fit[["sigma"]] )
std_resid_garch_norm <- fit3@fit[["residuals"]] / sqrt(fit3@fit[["sigma"]] )
std_resid_garch_t <- fit4@fit[["residuals"]] / sqrt(fit4@fit[["sigma"]] )

# Conduct LBQ test for each model (21 lags)
lbq_gjr_norm <- Box.test(std_resid_gjr_norm^2, lag = 21, type = "Ljung-Box")
lbq_gjr_t <- Box.test(std_resid_gjr_t^2, lag = 21, type = "Ljung-Box")
lbq_garch_norm <- Box.test(std_resid_garch_norm^2, lag = 21, type = "Ljung-Box")
lbq_garch_t <- Box.test(std_resid_garch_t^2, lag = 21, type = "Ljung-Box")

#test stat
lbq_gjr_norm$statistic
lbq_gjr_t$statistic
lbq_garch_norm$statistic
lbq_garch_t$statistic

# p value
lbq_gjr_norm$p.value
lbq_gjr_t$p.value
lbq_garch_norm$p.value
lbq_garch_t$p.value

# e) out of sample period, 1 step ahead forecasts ---- 
# Set the rolling window size
window_size <- 250

# Number of forecasts
n_forecasts <- length(data$log.return) - window_size

# Initialize vectors to store forecasts and squared errors
forecasts_gjr_norm <- numeric(n_forecasts)
forecasts_gjr_t <- numeric(n_forecasts)
forecasts_garch_norm <- numeric(n_forecasts)
forecasts_garch_t <- numeric(n_forecasts)

squared_errors_gjr_norm <- numeric(n_forecasts)
squared_errors_gjr_t <- numeric(n_forecasts)
squared_errors_garch_norm <- numeric(n_forecasts)
squared_errors_garch_t <- numeric(n_forecasts)

# Rolling window forecasting
for (i in 1:n_forecasts) {
  # Extract window data
  window_data <- data$log.return[(i):(i + window_size - 1)]
  
  # Fit models to window data
  gjr_fit_norm <- ugarchfit(spec = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1, 1)),
                                              mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                                              distribution.model = "norm"),
                            data = window_data)
  
  gjr_fit_t <- ugarchfit(spec = ugarchspec(variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1, 1)),
                                           mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                                           distribution.model = "std"),
                         data = window_data)
  
  garch_fit_norm <- ugarchfit(spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                                mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                                                distribution.model = "norm"),
                              data = window_data)
  
  garch_fit_t <- ugarchfit(spec = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                                             mean.model = list(armaOrder = c(0, 0), include.mean = TRUE),
                                             distribution.model = "std"),
                           data = window_data)
  
  # Produce 1-step ahead forecasts
  forecasts_gjr_norm[i] <- ugarchforecast(gjr_fit_norm, n.ahead = 1)@forecast$variance
  forecasts_gjr_t[i] <- ugarchforecast(gjr_fit_t, n.ahead = 1)@forecast$variance
  forecasts_garch_norm[i] <- ugarchforecast(garch_fit_norm, n.ahead = 1)@forecast$variance
  forecasts_garch_t[i] <- ugarchforecast(garch_fit_t, n.ahead = 1)@forecast$variance
  
  # Calculate squared errors
  squared_errors_gjr_norm[i] <- (returns[i + window_size] - 0)^2 - forecasts_gjr_norm[i]^2
  squared_errors_gjr_t[i] <- (returns[i + window_size] - 0)^2 - forecasts_gjr_t[i]^2
  squared_errors_garch_norm[i] <- (returns[i + window_size] - 0)^2 - forecasts_garch_norm[i]^2
  squared_errors_garch_t[i] <- (returns[i + window_size] - 0)^2 - forecasts_garch_t[i]^2
}

# Calculate RMSFE
rmsfe_gjr_norm <- sqrt(mean(squared_errors_gjr_norm^2))
rmsfe_gjr_t <- sqrt(mean(squared_errors_gjr_t^2))
rmsfe_garch_norm <- sqrt(mean(squared_errors_garch_norm^2))
rmsfe_garch_t <- sqrt(mean(squared_errors_garch_t^2))

# Print RMSFE
cat("RMSFE:\n")
cat("GJR-GARCH with Normal Innovations:", rmsfe_gjr_norm, "\n")
cat("GJR-GARCH with Student's t Innovations:", rmsfe_gjr_t, "\n")
cat("GARCH with Normal Innovations:", rmsfe_garch_norm, "\n")
cat("GARCH with Student's t Innovations:", rmsfe_garch_t, "\n")

# Apply DM test with GARCH(1, 1) model as benchmark
dm_test_gjr_norm <- dm.test(squared_errors_garch_norm, squared_errors_gjr_norm, power = 2, alternative = "less")
dm_test_gjr_t <- dm.test(squared_errors_garch_norm, squared_errors_gjr_t, power = 2, alternative = "less")
dm_test_garch_t <- dm.test(squared_errors_garch_norm, squared_errors_garch_t, power = 2, alternative = "less")

# Print DM test results
cat("\nDiebold-Mariano Test (GARCH(1, 1) as Benchmark):\n")
print(dm_test_gjr_norm)
print(dm_test_gjr_t)
print(dm_test_garch_t)



# claude  not classa

par(mfrow = c(1, 1))
plot(e, nicG1, type = "l", col = "blue", ylab = 'NIC', xlab = 'Shocks (e)', lwd = 2, lty = 1, main = "News Impact Curve")
lines(e, nicG2, type = "l", col = "red", lwd = 2, lty = 2)
lines(e, nicGJR1, type = "l", col = "green", lwd = 2, lty = 3)
lines(e, nicGJR2, type = "l", col = "purple", lwd = 2, lty = 4)
legend("top", legend = c("GARCH N(0,1)", "GARCH t_v", "GJR-GARCH N(0,1)", "GJR-GARCH t_v"), col = c("blue", "red", "green", "purple"), lty = 1:4, ncol = 1, lwd = 2)





