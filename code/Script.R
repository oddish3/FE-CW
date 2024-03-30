# ==============================================================================
# FInancaial Economitrics Coursework
# ==============================================================================
#
# Author: 10710007
# Version: 13-03-2024
#
# ==============================================================================

rm(list = ls())
# Packages
library(lubridate)
library(forecast)
library(rugarch)
# Data
setwd("/home/oddish3/Documents/R_folder/MSc/FE/FE-coursework/code")
data = read.csv("../data/group_11.csv")
source("fineco_fun.R")
source("../utils/latex-macro.R")
figures_path <- "../docs/figures/"


# Script
# ==============================================================================
data$date = as.Date(data$date)
# plot(data$log.return, type = "l", col = "darkgreen")
# abline(h = 0, v = 250, col = "red")
# dev.off()


# ------------------------------------------
#           Practical Exercise
# ------------------------------------------

# Use the first 250 observations as an in-sample (estimation) period and the last 250 observations as out of sample forecasting
in_sample = as.matrix(data[1:250, 2])
out_sample = as.matrix(data[(251:nrow(data)), 2])

# a) investigating statistical properties of the in-sample data ----
# i) descriptive stats
a1_results = dstats(in_sample)

# ii-v) moment tests
a2_results = test_moment(in_sample)

# Assemble data for moments and test statistics into a dataframe
moment_test_df = data.frame(
  # descriptive stats
  amu = a1_results[1,1],
  asigma = a1_results[2,1],
  askew = a1_results[3,1],
  akurt = a1_results[4,1],
  # t stats
  amut = a2_results[1,1],
  askewt = a2_results[2,1],
  akurtt = a2_results[3,1],
  ajbt = a2_results[4,1],
  # p vals
  amup = a2_results[1,2],
  askewp = a2_results[2,2],
  akurtp = a2_results[3,2],
  ajbp = a2_results[4,2]
)
# Appending the second section with its title
write_latex("../results/results.tex", moment_test_df, decimal_precision = 2, append = FALSE, section_title = "Moment Test and Descriptive Statistics")

# vi) lbq test
lbq1 = Box.test(in_sample, lag = 21, type = "Ljung-Box", fitdf = 0)
lbq2 = Box.test(in_sample^2, lag = 21, type = "Ljung-Box", fitdf = 0) # fitdf is the number of parameters estimated ???

lbq_df = data.frame(
  aistat = lbq1$statistic,
  aip = lbq1$p.value,
  aiistat = lbq2$statistic, 
  aiip = lbq2$p.value
)
# Writing the first section with its title
write_latex("../results/results.tex", lbq_df, decimal_precision = 2, append = TRUE, section_title = "Ljung-Box Test Results")

# vii) Plot SACF and SPACF
# Open a PNG device for ACF and PACF plots
png(paste0(figures_path, "PACF.png"))

# Setting the plotting area to accommodate two plots side by side
par(mfrow = c(1, 2))

# Generate ACF and PACF plots
Acf(in_sample, lag.max = 25)
Pacf(in_sample, lag.max = 25)

# Close the plotting device for ACF/PACF
dev.off()

# Reset par settings to default for subsequent plots
par(mfrow = c(1, 1))

# Now, generate and save another plot separately if needed
png(paste0(figures_path, "log_return_plot.png"))
plot(data$log.return, type = "l", col = "darkgreen")
dev.off()

# b) estimating the conditional variance ----
# Assume that the conditional mean of the return series is constant
#Use the et series to estimate the following conditional variance

# i) GARCH(1,1) with zt ~ N(0,1)
spec1 =  ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)), 
                    mean.model = list(armaOrder = c(0, 0), include.mean = FALSE) ,distribution.model = "norm")

fit1 = ugarchfit(spec = spec1, data = data$log.return)

estimates1 = fit1@fit$robust.matcoef[,1]  # This extracts the "Estimate" column
p_values1 = fit1@fit$robust.matcoef[,4]  # This extracts the "Pr(>|t|)" column


# GARCH(1,1) with zt ~ tv
spec2 = ugarchspec(variance.model = list(model = "sGARCH", garchOrder = c(1, 1)),
                   mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
                   distribution.model = "std")  # 'std' for Student's t-distribution
fit2 = ugarchfit(spec = spec2, data = data$log.return)

estimates2 = fit2@fit$robust.matcoef[,1]  # This extracts the "Estimate" column
p_values2 = fit2@fit$robust.matcoef[,4]  # This extracts the "Pr(>|t|)" column

# ii) GJR-GARCH (1, 1) with zt ~ N(0,1) 
spec3 = ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1,1)),
  mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
  distribution.model = "norm" # Standard normal distribution for innovations
)

fit3 = ugarchfit(spec = spec3, data = data$log.return)
estimates3 = fit3@fit$robust.matcoef[,1]  # This extracts the "Estimate" column
p_values3 = fit3@fit$robust.matcoef[,4]  # This extracts the "Pr(>|t|)" column

# GJR-MODEL (1, 1) with zt ~ tv
spec4 = ugarchspec(
  variance.model = list(model = "gjrGARCH", garchOrder = c(1, 1)),
  mean.model = list(armaOrder = c(0, 0), include.mean = FALSE),
  distribution.model = "std"  # 'std' for Student's t-distribution
)
fit4 = ugarchfit(spec = spec4, data = data$log.return)
estimates4 = fit4@fit$robust.matcoef[,1]  # This extracts the "Estimate" column
p_values4 = fit4@fit$robust.matcoef[,4] 

garch11_df = data.frame(
  bw = estimates1["omega"],
  ba = estimates1["alpha1"],
  bb = estimates1["beta1"],
  bpi = p_values1["omega"],
  bpii = p_values1["alpha1"],
  bpiii = p_values1["beta1"]
)

write_latex("../results/results.tex", garch11_df, append = TRUE, section_title = "GARCH(1,1) with zt ~ N(0,1)")

garch11_t_df = data.frame(
  bwi = estimates2["omega"],
  bai = estimates2["alpha1"],
  bbi = estimates2["beta1"],
  bvi = estimates2["shape"],
  bpti = p_values2["omega"],
  bptii = p_values2["alpha1"],
  bptiii = p_values2["beta1"],
  bptiv = p_values2["shape"]
)
write_latex("../results/results.tex", garch11_t_df, append = TRUE, section_title = "GARCH(1,1) with zt ~ tv")

garch_gjr_df = data.frame(
  bwii = estimates3["omega"],
  baii = estimates3["alpha1"],
  bbii = estimates3["beta1"],
  bgii = estimates3["gamma1"],
  bpgi = p_values3["omega"],
  bpgii = p_values3["alpha1"],
  bpgiii = p_values3["beta1"],
  bpgiv = p_values3["gamma1"]
)

write_latex("../results/results.tex", garch_gjr_df, append = TRUE, section_title = "GJR-GARCH(1,1) with zt ~ N(0,1)")

garch_gjr_t_df = data.frame(
  bwiii = estimates4["omega"],
  baiii = estimates4["alpha1"],
  bbiii = estimates4["beta1"],
  bgiii = estimates4["gamma1"],
  bviii = estimates4["shape"],
  bpgtp = p_values4["omega"],
  bpgtpi = p_values4["alpha1"],
  bpgtpii = p_values4["beta1"],
  bpgtpiii = p_values4["gamma1"],
  bpgtpiv = p_values4["shape"]
)

write_latex("../results/results.tex", garch_gjr_t_df, append = TRUE, section_title = "GJR-GARCH(1,1) with zt ~ tv")

# c)  plotting NIC ----
# For GARCH(1,1) with zt ~ N(0,1)
w1 = estimates1["omega"]
a1 = estimates1["alpha1"]
b1 = estimates1["beta1"]

# For GARCH(1,1) with zt ~ tv
w2 = estimates2["omega"]
a2 = estimates2["alpha1"]
b2 = estimates2["beta1"]
v2 = estimates2["shape"]

# For GJR-GARCH(1,1) with zt ~ N(0,1)
w3 = estimates3["omega"]
a3 = estimates3["alpha1"]
b3 = estimates3["beta1"]
g3 = estimates3["gamma1"]

# For GJR-GARCH(1,1) with zt ~ tv
w4 = estimates4["omega"]
a4 = estimates4["alpha1"]
b4 = estimates4["beta1"]
g4 = estimates4["gamma1"]
v4 = estimates4["shape"]

# Unconditional variance from the first GARCH(1,1) model
ve = w1 / (1 - a1 - b1)

# NIC from tutorial ----
T = 500
e = seq(-5, 5, length.out = T)  # Grid of shocks epsilon
nicG1 = nicG2 = nicGJR1 = nicGJR2 = rep(0, T)  # Initialize NIC for each model

# Calculate NIC for each model
for (t in 1:T) {
  nicG1[t] = w1 + b1 * ve + a1 * e[t]^2  # GARCH(1,1) with zt ~ N(0,1)
  nicG2[t] = w2 + b2 * ve + a2 * e[t]^2  # GARCH(1,1) with zt ~ tv

  if (e[t] > 0) { 
    nicGJR1[t] = w3 + b3 * ve + a3 * e[t]^2 # GJR-GARCH(1,1) with zt ~ N(0,1)
    nicGJR2[t] = w4 + b4 * ve + a4 * e[t]^2 # GJR-GARCH(1,1) with zt ~ tv
  }
  else{
    nicGJR1[t] = w3 + b3 * ve + (a3 + g3) * e[t]^2 # GJR-GARCH(1,1) with zt ~ N(0,1)
    nicGJR2[t] = w4 + b4 * ve + (a4 + g4) * e[t]^2 # GJR-GARCH(1,1) with zt ~ tv
  }
}

# Plot NIC for all four models
par(mfrow = c(1, 1))
png(paste0(figures_path, "NIC.png"))
plot(e, nicG1, type = "l", col = "blue", ylab = 'NIC', xlab = 'Shocks (e)', lwd = 2, lty = 1, main = "News Impact Curve")
lines(e, nicG2, type = "l", col = "red", lwd = 2, lty = 2)
lines(e, nicGJR1, type = "l", col = "green", lwd = 2, lty = 3)
lines(e, nicGJR2, type = "l", col = "purple", lwd = 2, lty = 4)
legend("top", legend = c("GARCH N(0,1)", "GARCH t_v", "GJR-GARCH N(0,1)", "GJR-GARCH t_v"), col = c("blue", "red", "green", "purple"), lty = 1:4, ncol = 1, lwd = 2)
dev.off()

# d) best model by analysing standarised residuals ----

# Extract standardized residuals from each model
# Conduct LBQ test for each model (21 lags) residuals
# Conduct LBQ test for each model (21 lags) squared residuals

# For GARCH(1,1) with zt ~ N(0,1)
std_resid_norm = fit1@fit[["residuals"]] / sqrt(fit1@fit[["sigma"]] )
lbq_z_norm = Box.test(std_resid_norm, lag = 21, type = "Ljung-Box")
lbq_z_norm2 = Box.test(std_resid_norm^2, lag = 21, type = "Ljung-Box")

# For GARCH(1,1) with zt ~ tv
std_resid_t = fit2@fit[["residuals"]] / sqrt(fit2@fit[["sigma"]] )
lbq_z_t = Box.test(std_resid_t, lag = 21, type = "Ljung-Box")
lbq_z_t2 = Box.test(std_resid_t^2, lag = 21, type = "Ljung-Box")

# For GJR-GARCH(1,1) with zt ~ N(0,1)
std_resid_gjr_norm = fit3@fit[["residuals"]] / sqrt(fit3@fit[["sigma"]] )
lbq_z_gjr_norm = Box.test(std_resid_gjr_norm, lag = 21, type = "Ljung-Box")
lbq_z_gjr_norm2 = Box.test(std_resid_gjr_norm^2, lag = 21, type = "Ljung-Box")

# For GJR-GARCH(1,1) with zt ~ tv
std_resid_gjr_t = fit4@fit[["residuals"]] / sqrt(fit4@fit[["sigma"]] )
lbq_z_gjr_t = Box.test(std_resid_gjr_t, lag = 21, type = "Ljung-Box")
lbq_z_gjr_t2 = Box.test(std_resid_gjr_t^2, lag = 21, type = "Ljung-Box")


df_test = data.frame(
  # garch N(0,1)
  zone = lbq_z_norm$statistic,
  pone = lbq_z_norm$p.value,
  zfive = lbq_z_norm2$statistic,
  pfive = lbq_z_norm2$p.value,
  # garch t
  ztwo = lbq_z_t$statistic,
  ptwo = lbq_z_t$p.value,
  zsix = lbq_z_t2$statistic,
  psix = lbq_z_t2$p.value,
  # gjr N(0,1)
  zthree = lbq_z_gjr_norm$statistic,
  pthree = lbq_z_gjr_norm$p.value,
  zseven = lbq_z_gjr_norm2$statistic,
  pseven = lbq_z_gjr_norm2$p.value,
  # gjr t
  zfour = lbq_z_gjr_t$statistic,
  pfour = lbq_z_gjr_t$p.value,
  zeight = lbq_z_gjr_t2$statistic,
  peight = lbq_z_gjr_t2$p.value
)

write_latex("../results/results.tex", df_test, append = TRUE, section_title = "residuals and squared lbq")

# e) 1 step 2ahead forecasting the conditional variance ----

H = 250
T = length(data$log.return) - H
f1 = f2 = f3 = f4 = matrix(0, H, 1)  # Initialize forecast for each model

for (i in 1:H) {
  window = data$log.return[i:(T+i-1)]
  
  fit.g11 = ugarchfit(spec = spec1, data = window, solver = 'hybrid')
  fit.tg11 = ugarchfit(spec = spec2, data = window, solver = 'hybrid')
  fit.gj11 = ugarchfit(spec = spec3, data = window, solver = 'hybrid')
  fit.tgj11 = ugarchfit(spec = spec4, data = window, solver = 'hybrid')
  
  # forecast
  xx = ugarchforecast(fit.g11, data = window, n.ahead = 1)
  f1[i] = xx@forecast$sigmaFor
  
  xx = ugarchforecast(fit.tg11, data = window, n.ahead = 1)
  f2[i] = xx@forecast$sigmaFor
  
  xx = ugarchforecast(fit.gj11, data = window, n.ahead = 1)
  f3[i] = xx@forecast$sigmaFor
  
  xx = ugarchforecast(fit.tgj11, data = window, n.ahead = 1)
  f4[i] = xx@forecast$sigmaFor
  
  print(i)
}

# Forecast errors
e1 = f1 - data$log.return[(T+1):(T+H)]^2
e2 = f2 - data$log.return[(T+1):(T+H)]^2
e3 = f3 - data$log.return[(T+1):(T+H)]^2
e4 = f4 - data$log.return[(T+1):(T+H)]^2

# RMSFE
rmsfe = function(e) {
  sse = sum(e^2) / length(e)
  r = sqrt(sse)
  return(r)
}

# dm test
dm_test_12 = dm.test(e1, e2, alternative = "less", h = 1, power = 2, varestimator = "acf")
dm_test_13 = dm.test(e1, e3, alternative = "less", h = 1, power = 2, varestimator = "acf")
dm_test_14 = dm.test(e1, e4, alternative = "less", h = 1, power = 2, varestimator = "acf")

# Extract DM test statistics and p-values
dm_stat_12 = dm_test_12$statistic
p_value_12 = dm_test_12$p.value

dm_stat_13 = dm_test_13$statistic
p_value_13 = dm_test_13$p.value

dm_stat_14 = dm_test_14$statistic
p_value_14 = dm_test_14$p.value


rmsfe_dm_df = data.frame(
  #rmsfe
  rmsfei = rmsfe(e1),
  rmsfeii = rmsfe(e2),
  rmsfeiii = rmsfe(e3),
  rmsfeiv = rmsfe(e4),
  #dm test
  dm = dm_stat_12,
  dmpii = p_value_12,
  dmi = dm_stat_13,
  dmpiv = p_value_13,
  dmii = dm_stat_14,
  dmpv = p_value_14
  
)
write_latex("../results/results.tex", rmsfe_dm_df, append = TRUE, section_title = "root mean square forecast error and dm test")



# End of Script
# ==============================================================================
