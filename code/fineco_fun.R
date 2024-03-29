lr = function(p){
  r = log(p[-1]) - log(head(p,-1))
  return(r)
}

dstats = function(x){
  T = nrow(x)
  'Sample mean'
  mu = colSums(x)/T
  'Sample standard deviation'
  sigma2 = colSums((x-mu)^2)/T
  sigma = sqrt(sigma2)
  'Sample skewness'
  skew = colSums((x-mu)^3/(sigma^3))/T
  'Sample kurtosis'
  kurt = colSums((x-mu)^4/(sigma^4))/T
  d = rbind(mu,sigma,skew,kurt)
  return(d)
}


test_moment = function(x){
  x = as.matrix(x)
  T = nrow(x)
  moments = dstats(x)
  'Testing the mean against zero'
  t_mean = sqrt(T)*moments[1]/moments[2]
  pv_mean =2*(1- pnorm(abs(t_mean)))
  'Testing skewness against zero'
  t_skew = moments[3]/sqrt(6/T)
  pv_skew =2*(1- pnorm(abs(t_skew)))
  'Testing kurtosis against three'
  t_kurt = (moments[4]-3)/sqrt(24/T)
  pv_kurt =2*(1- pnorm(abs(t_kurt)))
  'Jarque-Bera test for normality'
  t_JB = (T/6)*(moments[3]**2 + ((moments[4]-3)**2)/4)
  pv_JB =(1- pchisq(t_JB,2))
  tstat = c(t_mean, t_skew, t_kurt, t_JB)
  pv = c(pv_mean, pv_skew, pv_kurt, pv_JB)
  return (cbind(tstat,pv))
}

# JB test function for matrices
JBtest = function(x){
  T = nrow(x)
  moments = dstats(x)
  #Jarque-Bera test for normality
  t_JB = (T/6)*(moments[3,]**2 + ((moments[4,]-3)**2)/4)
  pv_JB =(1- pchisq(t_JB,2))
  return (rbind(t_JB,pv_JB))
}

# ARCH-LM test by Engle: short function
archlm = function(e,lag){
  archlm = ar.ols(e^2, aic = F, order.max = lag, intercept = T)
  # save residuals
  eps = archlm$resid
  eps = na.omit(eps)
  # Conduct an LM test
  R2 = 1-sum(eps^2)/sum((e^2-mean(e^2))^2)
  LM = length(e)*R2
  lm_pvalue = 1-pchisq(LM,lag)
  return(lm_pvalue)
}