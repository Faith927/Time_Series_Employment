rm(list=ls())
library(TSA)
library(tseries)
library(forecast)
library(lmtest)
library(dplyr)

residual.analysis <- function(model, std = TRUE,start = 2, class = c("ARIMA","GARCH","ARMA-GARCH", "garch", "fGARCH")[1]){
  library(TSA)
  library(FitAR)
  if (class == "ARIMA"){
    if (std == TRUE){
      res.model = rstandard(model)
    }else{
      res.model = residuals(model)
    }
  }else if (class == "GARCH"){
    res.model = model$residuals[start:model$n.used]
  }else if (class == "garch"){
    res.model = model$residuals[start:model$n.used]  
  }else if (class == "ARMA-GARCH"){
    res.model = model@fit$residuals
  }else if (class == "fGARCH"){
    res.model = model@residuals
  }else {
    stop("The argument 'class' must be either 'ARIMA' or 'GARCH' ")
  }
  
  par(mfrow=c(3,2))
  plot(res.model,type='o',ylab='Standardised residuals', main="Time series plot of standardised residuals")
  abline(h=0)
  hist(res.model,main="Histogram of standardised residuals")
  qqnorm(res.model,main="QQ plot of standardised residuals")
  qqline(res.model, col = 2)
  acf(res.model,main="ACF of standardised residuals")
  print(shapiro.test(res.model))
  k=0
  LBQPlot(res.model, lag.max = 30, StartLag = k + 1, k = 0, SquaredQ = FALSE)
  par(mfrow=c(1,1))
}

show_ACF_PACF <- function(series,maxlag,name){
  par(mfrow=c(1,2))
  seasonal_acf(series,lag.max=maxlag, main=paste("ACF plot of",name))
  seasonal_pacf(series, lag.max=maxlag, main=paste("PACF plot of",name))
  par(mfrow=c(1,1))
}


sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}


# Normalize function
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}


helper <- function(class = c("acf", "pacf"), ...) {
  
  # Capture additional arguments
  params <- match.call(expand.dots = TRUE)
  params <- as.list(params)[-1]
  
  # Calculate ACF/PACF values
  if (class == "acf") {
    acf_values <- do.call(acf, c(params, list(plot = FALSE)))
  } else if (class == "pacf") {
    acf_values <- do.call(pacf, c(params, list(plot = FALSE)))
  }
  
  # Extract values and lags
  acf_data <- data.frame(
    Lag = as.numeric(acf_values$lag),  
    ACF = as.numeric(acf_values$acf)   
  )
  
  # Identify seasonal lags to be highlighted
  seasonal_lags <- acf_data$Lag %% 1 == 0
  
  # Plot ACF/PACF values
  if (class == "acf") {
    do.call(acf, c(params, list(plot = TRUE)))
  } else if (class == "pacf") {
    do.call(pacf, c(params, list(plot = TRUE)))
  }
  
  # Add colored segments for seasonal lags
  for (i in which(seasonal_lags)) {
    segments(x0 = acf_data$Lag[i], y0 = 0, x1 = acf_data$Lag[i], y1 = acf_data$ACF[i], col = "red")
  }
}


# seasonal_acf ------------------------------------------------------------

seasonal_acf <- function(...) {
  helper(class = "acf", ...)
}


# seasonal_pacf -----------------------------------------------------------

seasonal_pacf <- function(...) {
  helper(class = "pacf", ...)
}

show_ACF_PACF <- function(series,maxlag,name){
  par(mfrow=c(1,2))
  seasonal_acf(series,lag.max=maxlag, main=paste("ACF plot of",name))
  seasonal_pacf(series, lag.max=maxlag, main=paste("PACF plot of",name))
  par(mfrow=c(1,1))
}


# read in csv file
employment <- read.csv("employment.csv", header=TRUE)
employment
class(employment) 
head(employment)


# Conversion to Time Series
employment_TS <- ts(employment$employed, start = c(1978, 2), frequency = 12)
#Subset the series ot 2010-2024 only
subset = window(employment_TS,2010,)
subset

class(subset)
summary(subset)
x11()
plot(subset,type='l',ylab='Employment', main = " Time series plot of Employment by 1000")
summary(subset)

# acf and pacf plots
show_ACF_PACF(subset,maxlag = 48,name = "the Employment series.")


# adf test (weak indicator of non-stationary)
adf.test(subset, alternative = c("stationary"))
# pp test
pp.test(subset)
# kpss test
kpss.test(subset)

# # MCleod test
# McLeod.Li.test(y=employment_TS,main="McLeod-Li Test Statistics for Monthly Employment Numbers")

# qq plot 
qqnorm(subset)
qqline(subset, col = 2)
shapiro.test(subset)

# Normalization 
# A small positive constant
log_employment <- log(subset + .001)
min_value <- min(subset)

# qq plot
qqnorm(log_employment)
qqline(log_employment, col = 2)
shapiro.test(log_employment)


normalized_ts <- normalize(subset)

# qq plot 
qqnorm(normalized_ts)
qqline(normalized_ts, col = 2)
shapiro.test(normalized_ts)


# Box Cox Transformation
BC <- BoxCox.ar(subset, lambda = seq(0.5, 1, 0.01))
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda # lambda = 0.5 corresponds to square root transformation?

sqrt_data <- sqrt(subset)

plot(sqrt_data,type='l',ylab= 'Employment Numbers', main = 'Time series plot of BC transformed 
     monthly average employment numbers.')

# acf and pacf plots
par(mfrow=c(1,2))
show_ACF_PACF(sqrt_data,maxlag = 48,name = "the Employment series.")
par(mfrow=c(1,1))

# differencing

diff.employment = diff(sqrt_data)

par(mfrow=c(1,2))
show_ACF_PACF(diff.employment,maxlag=48,name="Squared-rooted subset series")
par(mfrow=c(1,1))

# adf, pp, kpss test
adf.test(diff.employment)
pp.test(diff.employment) 
kpss.test(diff.employment)

#SEASONAL ARIMA
sqrt_data
#Start with a baseline model with D=1 
m1 = Arima(sqrt_data,order=c(0,0,0),seasonal=list(order=c(0,1,0), period=12))
res.m1 = residuals(m1)
par(mfrow=c(1,1))
plot(res.m1,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
show_ACF_PACF(res.m1, maxlag = 78, name="Residuals of SARIMA(0,0,0)x(0,1,0)_12")

#From the ACF and PACF => P=0 and Q=2
m2 = Arima(sqrt_data,order=c(0,0,0),seasonal=list(order=c(0,1,1), period=12),method="CSS")
res.m2 = residuals(m2)
par(mfrow=c(1,1))
plot(res.m2,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
show_ACF_PACF(res.m2, maxlag = 78,name="Residuals of SARIMA(0,0,0)x(0,1,1)_12")
adf.test(res.m2) #Residuals are stationary

#Add differencing
m3 = Arima(sqrt_data,order=c(0,1,0),seasonal=list(order=c(0,1,1), period=12))
res.m3 = residuals(m3);  
par(mfrow=c(1,1))
plot(res.m3,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
show_ACF_PACF(res.m3, maxlag = 48,name="Residuals of SARIMA(0,1,0)x(0,1,1)_12")
#No trend or seasonality left, move forward
adf.test(res.m3)

#From ACF and PACF of m3
#p = 3 and q = 2
#Determined by significant bars between 0 and 1
m4 = Arima(sqrt_data,order=c(3,1,2),
           seasonal=list(order=c(0,1,1), period=12))
res.m4 = residuals(m4);  
par(mfrow=c(1,1))
plot(res.m4,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
show_ACF_PACF(res.m4, maxlag = 48,name="Residuals of SARIMA(3,1,3)x(0,1,2)_12")
# SARIMA(3,1,2)x(0,1,1)_12

eacf(res.m3) 
# The tentative models are specified as 
# SARIMA(0,1,2)x(0,1,1)_12
# SARIMA(0,1,3)x(0,1,1)_12
# SARIMA(1,1,2)x(0,1,1)_12
# SARIMA(1,1,3)x(0,1,1)_12

par(mfrow=c(1,1))
bic_table = armasubsets(y=res.m3,nar=15,nma=15,y.name='p',ar.method='ols')
plot(bic_table)
# SARIMA(0,1,14)x(0,1,1)_12


#Final sets of possible models
# SARIMA(3,1,2)x(0,1,1)_12
# SARIMA(0,1,2)x(0,1,1)_12
# SARIMA(0,1,3)x(0,1,1)_12
# SARIMA(1,1,2)x(0,1,1)_12
# SARIMA(1,1,3)x(0,1,1)_12
# SARIMA(0,1,14)x(0,1,1)_12

#Define a function to fit models quickly
fit_SARIMA <- function(series,orders) {
  model_css = Arima(series,order=orders,
                    seasonal=list(order=c(0,1,1),period=12),
                    method='CSS')
  coef_css = coeftest(model_css)
  print("CSS Method")
  print(coef_css)
  
  model_ml = Arima(series,order=orders,
                   seasonal=list(order=c(0,1,1),period=12),
                   method='ML')
  coef_ml = coeftest(model_ml)
  print("ML Method")
  print(coef_ml)
  return(list(css=model_css,ml=model_ml))
}
# SARIMA(3,1,2)x(0,1,1)_12
m312 = fit_SARIMA(sqrt_data,orders=c(3,1,2))
residual.analysis(model = m312$css) #Did not capture the info well
# SARIMA(0,1,2)x(0,1,1)_12
m012 = fit_SARIMA(sqrt_data,orders=c(0,1,2))
residual.analysis(model = m012$css) #Did not capture the info well
# SARIMA(0,1,3)x(0,1,1)_12
m013 = fit_SARIMA(sqrt_data,orders=c(0,1,3))
residual.analysis(model = m013$css) #Did not capture the info well
# SARIMA(1,1,2)x(0,1,1)_12
m112 = fit_SARIMA(sqrt_data,orders=c(1,1,2))
residual.analysis(model = m112$css) #Did not capture the info well
# SARIMA(1,1,3)x(0,1,1)_12
m113 = fit_SARIMA(sqrt_data,orders=c(1,1,3))
residual.analysis(model = m113$css) #Did not capture the info well
# SARIMA(0,1,14)x(0,1,1)_12
m0114 = fit_SARIMA(sqrt_data,orders=c(0,1,14))
residual.analysis(model = m0114$css) #Captured info pretty well


sc.AIC = AIC(m312$ml,m012$ml,m013$ml,m112$ml,m113$ml,m0114$ml)
sc.BIC = BIC(m312$ml,m012$ml,m013$ml,m112$ml,m113$ml,m0114$ml)
sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "bic")
#The big model captured information well but have the worst aic and bic scores

# error metrics for first 6 models:
#ERROR METRICS
m312css <- accuracy(m312$css)[1:7]
m012css <- accuracy(m012$css)[1:7]
m013css <- accuracy(m013$css)[1:7]
m112css <- accuracy(m112$css)[1:7]
m113css <- accuracy(m113$css)[1:7]
m0114css <- accuracy(m0114$css)[1:7]

df.Smodels <- data.frame(
  rbind(m312css, m012css, m013css, m112css, m113css, m0114css)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("SARIMA(3,1,2)x(0,1,1)_12", "SARIMA(0,1,2)x(0,1,1)_12", "SARIMA(0,1,3)x(0,1,1)_12", 
                          "SARIMA(1,1,2)x(0,1,1)_12", "SARIMA(1,1,3)x(0,1,1)_12", "SARIMA(0,1,14)x(0,1,1)_12"
)
round(df.Smodels,  digits = 3)


#Residuals of best models from BIC still have
#significant autocorrelation
residual.analysis(model = m012$css)
#Best models from AIC captured the information well
#and nothing valuable is in their residuals

m1114 = fit_SARIMA(sqrt_data,orders=c(1,1,14)) #About half of the vars are significant
residual.analysis(m1114$ml) #Capture info pretty well

#Try overfitted model
m0115 = fit_SARIMA(sqrt_data,orders=c(0,1,15)) #Most vars are significant in CSS
residual.analysis(m0115$ml) #Capture info very well


#ERROR METRICS
m312css <- accuracy(m312$css)[1:7]
m012css <- accuracy(m012$css)[1:7]
m013css <- accuracy(m013$css)[1:7]
m112css <- accuracy(m112$css)[1:7]
m113css <- accuracy(m113$css)[1:7]
m0114css <- accuracy(m0114$css)[1:7]
m0115css <- accuracy(m0115$css)[1:7]
m1114css <- accuracy(m1114$css)[1:7]

df.Smodels <- data.frame(
  rbind(m312css, m012css, m013css, m112css, m113css, m0114css, m0115css,m1114css)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("SARIMA(3,1,2)x(0,1,1)_12", "SARIMA(0,1,2)x(0,1,1)_12", "SARIMA(0,1,3)x(0,1,1)_12", 
                          "SARIMA(1,1,2)x(0,1,1)_12", "SARIMA(1,1,3)x(0,1,1)_12", "SARIMA(0,1,14)x(0,1,1)_12",
                          "SARIMA(0,1,15)x(0,1,1)_12","SARIMA(1,1,14)x(0,1,1)_12")
round(df.Smodels,  digits = 3)


sc.AIC = AIC(m312$ml,m012$ml,m013$ml,m112$ml,m113$ml,m0114$ml,m0115$ml,m1114$ml)
sc.BIC = BIC(m312$ml,m012$ml,m013$ml,m112$ml,m113$ml,m0114$ml,m0115$ml,m1114$ml)
sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "bic")
#Big models capture information well but have the worst aic and bic scores

#AIC picks m316 over 315

#FORECASTING


FC_CSS = Arima(subset,order=c(0,1,15),seasonal=list(order=c(0,1,1), period=12), 
              lambda = 0.5, method = "CSS")
forecastCSS = forecast(FC_CSS,lambda = 0.5,  h = 10)
forecastCSS
plot(forecastCSS)

FC_ML = Arima(subset,order=c(0,1,15),seasonal=list(order=c(0,1,1), period=12), 
               lambda = 0.5, method = "ML")
forecastML = forecast(FC_ML,lambda = 0.5,  h = 10)
forecastML
plot(forecastML)

