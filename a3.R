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

sort.score <- function(x, score = c("bic", "aic")){
  if (score == "aic"){
    x[with(x, order(AIC)),]
  } else if (score == "bic") {
    x[with(x, order(BIC)),]
  } else {
    warning('score = "x" only accepts valid arguments ("aic","bic")')
  }
}

# read in csv file
employment <- read.csv("employment.csv", header=TRUE)
employment
class(employment) 
head(employment)


# Conversion to Time Series
employment_TS <- ts(employment$employed, start = c(1978, 2), frequency = 12)
class(employment_TS)
summary(employment_TS)
plot(employment_TS,type='o',ylab='Employment', main = " Time series plot of Employment by 1000")
summary(employment_TS)

# Calculate the Z-scores
z_scores <- scale(employment_TS)

# Identify outliers (e.g., Z-score > 3 or < -3)
outliers_z <- which(abs(z_scores) > 3)
outliers_z

# adf test (weak indicator of non-stationary)
adf.test(employment_TS, alternative = c("stationary")) # Use default value of k
# pp test
pp.test(employment_TS)
# kpss test
kpss.test(employment_TS)

# MCleod test
McLeod.Li.test(y=employment_TS,main="McLeod-Li Test Statistics for Daily Google Returns")

# acf and pacf plots
par(mfrow=c(1,2))
acf(employment_TS,lag.max = 98, main ="ACF plot for the Global Temperature Anomalies series.")
pacf(employment_TS, main ="PACF plot for the Global Temperature Anomalies series.")
par(mfrow=c(1,1))


# qq plot 
qqnorm(employment_TS)
qqline(employment_TS, col = 2)
shapiro.test(employment_TS)


# Box Cox Transformation
BC <- BoxCox.ar(employment_TS, lambda = seq(0.5, 1, 0.01))
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda # lambda = 0.5 corresponds to square root transformation?

sqrt_data <- sqrt(employment_TS)

plot(sqrt_data,type='o',ylab='Time series plot of BC transformed 
     yearly average employment numbers.')

# acf and pacf plots
par(mfrow=c(1,2))
acf(sqrt_data, main='ACF plot of employment series.')
pacf(sqrt_data, main='PACF plot of employment series.')
par(mfrow=c(1,1))

# differencing
diff.employment = diff(sqrt_data)
par(mfrow=c(1,1))
plot(diff.employment,type='o',ylab='Average employment numbers', main = "Time series plot of the first differenced
     yearly average employment numbers.")
abline(h=0)


# strong correlation and arch affect
McLeod.Li.test(y=diff.employment, main = "McLeod-Li Test Statistics for Monthly Employment Numbers")
qqnorm(diff.employment)
# adf, pp, kpss test, mcleod
adf.test(diff.employment) # stationary
pp.test(diff.employment) # stationary
kpss.test(diff.employment) # non stationary
McLeod.Li.test(y=diff.employment,main="McLeod-Li Test Statistics for Daily Google Returns")


hist(diff.employment, 
     main = "Histogram of diff.employment",
     xlab = "Difference in Employment",
     ylab = "Frequency",
     col = "skyblue",
     border = "white")

# qq plot 
qqnorm(diff.employment)
qqline(diff.employment, col = 2)
shapiro.test(diff.employment)

# pacf and acf plots
par(mfrow=c(1,2))
acf(diff.employment, main='ACF plot of unemployment series.')
pacf(diff.employment, main='PACF plot of unemployment series.')
par(mfrow=c(1,1))


# qq plot 
qqnorm(diff.employment)
qqline(diff.employment, col = 2)
shapiro.test(diff.employment)



# Normalization 
# A small positive constant
log_diff_employment <- log(diff.employment + 2.7)
min_value <- min(diff.employment)
min_value

# qq plot
qqnorm(log_diff_employment)
qqline(log_diff_employment, col = 2)
shapiro.test(log_diff_employment)



# Normalize function
normalize <- function(x) {
  return((x - min(x)) / (max(x) - min(x)))
}

normalized_ts <- normalize(diff.employment)

# qq plot 
qqnorm(normalized_ts)
qqline(normalized_ts, col = 2)
shapiro.test(normalized_ts)

# need arma model to capture trend of the series
# significatn lags in acf and pacf and mcleod significant so then do arma + garch
# ARCH

# module 9, increasing trend with a variation that get higher as time passes
# indicates a higher variability with higher levels of employment
# do log
diff.employment = diff(log(sqrt_data))
plot(diff.employment)
# strong correlation and arch affect
McLeod.Li.test(y=diff.employment, main = "McLeod-Li Test Statistics for Monthly Employment Numbers")
qqnorm(diff.employment)
# adf, pp, kpss test, mcleod
adf.test(diff.employment) # stationary
pp.test(diff.employment) # stationary
kpss.test(diff.employment) # non stationary
McLeod.Li.test(y=diff.employment,main="McLeod-Li Test Statistics for Monthly Employment Numbers")



par(mfrow=c(1,2))
acf(diff.employment, main="The ACF plot of returns series for Monthly Employment Numbers")
pacf(diff.employment, main="The PACF plot of returns series Monthly Employment Numbers")
par(mfrow=c(1,1))
# a lot of significant lags...
# ARIMA(7,1, 8), ARIMA(8, 1, 8), 

eacf(diff.employment)
# {ARMA(4,1,8), ARMA(4,1,9), ARMA(5,1,8), ARMA(5,1,8)}
# chose because more parameters to choose from and 5,6 is 11 anyway


res = armasubsets(y=diff.employment,nar=14,nma=14,y.name='p',ar.method='ols')
plot(res)

# {ARIMA(12, 1, 12), ARIMA(12, 1, 10), ARIMA (12, 1, 6)}

# Overall models:
# {ARIMA(7,1, 8), ARIMA(8, 1, 8), ARMA(4,1,8), ARMA(4,1,9), ARMA(5,1,8), 
# ARMA(5,1,9), ARIMA(12, 1, 12), ARIMA(12, 1, 10), ARIMA (12, 1, 6)}


# ARIMA(7,1,8)
model_718_css = Arima(diff.employment,order=c(7,1,8),method='CSS')
coeftest(model_718_css)
residual.analysis(model = model_718_css)

model_718_ml = Arima(diff.employment,order=c(7,1,8),method='ML')
coeftest(model_718_ml)
residual.analysis(model = model_718_ml)

model_718_ml_css = Arima(diff.employment,order=c(7,1,8),method='CSS-ML')
coeftest(model_718_ml_css)
residual.analysis(model = model_718_ml_css)



# ARIMA(8,1,8)
model_818_css = Arima(diff.employment,order=c(8,1,8),method='CSS')
coeftest(model_818_css)
residual.analysis(model = model_818_css)

model_818_ml = Arima(diff.employment,order=c(8,1,8),method='ML')
coeftest(model_818_ml)
residual.analysis(model = model_818_ml)

model_818_ml_css = Arima(diff.employment,order=c(8,1,8),method='CSS-ML')
coeftest(model_818_ml_css)
residual.analysis(model = model_818_ml_css)


# ARIMA(4,1,8)
model_418_css = Arima(diff.employment,order=c(4,1,8),method='CSS')
coeftest(model_418_css)
residual.analysis(model = model_418_css)

model_418_ml = Arima(diff.employment,order=c(4,1,8),method='ML')
coeftest(model_418_ml)
residual.analysis(model = model_418_ml)

model_418_ml_css = Arima(diff.employment,order=c(4,1,8),method='CSS-ML')
coeftest(model_418_ml_css)
residual.analysis(model = model_418_ml_css)

# ARIMA(4,1,9)
model_419_css = Arima(diff.employment,order=c(4,1,9),method='CSS')
coeftest(model_419_css)
residual.analysis(model = model_419_css)

model_419_ml = Arima(diff.employment,order=c(4,1,9),method='ML')
coeftest(model_419_ml)
residual.analysis(model = model_419_ml)

model_419_ml_css = Arima(diff.employment,order=c(4,1,9),method='CSS-ML')
coeftest(model_419_ml_css)
residual.analysis(model = model_419_ml_css)

# ARIMA(5,1,8)
model_518_css = Arima(diff.employment,order=c(5,1,8),method='CSS')
coeftest(model_518_css)
residual.analysis(model = model_518_css)

model_518_ml = Arima(diff.employment,order=c(5,1,8),method='ML')
coeftest(model_518_ml)
residual.analysis(model = model_518_ml)
# Doesn't work
# model_518_ml_css = Arima(diff.employment,order=c(5,1,8),method='CSS-ML')
# coeftest(model_518_ml_css)
# residual.analysis(model = model_518_ml_css)

# ARIMA(5,1,9)
model_519_css = Arima(diff.employment,order=c(5,1,9),method='CSS')
coeftest(model_519_css)
residual.analysis(model = model_519_css)

model_519_ml = Arima(diff.employment,order=c(5,1,9),method='ML')
coeftest(model_519_ml)
residual.analysis(model = model_519_ml)

model_519_ml_css = Arima(diff.employment,order=c(5,1,9),method='CSS-ML')
coeftest(model_519_ml_css)
residual.analysis(model = model_519_ml_css)


# ARIMA(12,1,12)
model_12112_css = Arima(diff.employment,order=c(12,1,12),method='CSS')
coeftest(model_12112_css)
residual.analysis(model = model_12112_css)

model_12112_ml = Arima(diff.employment,order=c(12,1,12),method='ML')
coeftest(model_12112_ml)
residual.analysis(model = model_12112_ml)
# doesn't work
# model_12112_ml_css = Arima(diff.employment,order=c(12,1,12),method='CSS-ML')
# coeftest(model_12112_ml_css)
# residual.analysis(model = model_12112_ml_css)

# ARIMA(12,1,10)
model_12110_css = Arima(diff.employment,order=c(12,1,10),method='CSS')
coeftest(model_12110_css)
residual.analysis(model = model_12110_css)

model_12110_ml = Arima(diff.employment,order=c(12,1,10),method='ML')
coeftest(model_12110_ml)
residual.analysis(model = model_12110_ml)

model_12110_ml_css = Arima(diff.employment,order=c(12,1,10),method='CSS-ML')
coeftest(model_12110_ml_css)
residual.analysis(model = model_12110_ml_css)

# ARIMA(12,1,6)
model_1216_css = Arima(diff.employment,order=c(12,1,6),method='CSS')
coeftest(model_1216_css)
residual.analysis(model = model_1216_css)

model_1216_ml = Arima(diff.employment,order=c(12,1,6),method='ML')
coeftest(model_1216_ml)
residual.analysis(model = model_1216_ml)

model_1216_ml_css = Arima(diff.employment,order=c(12,1,6),method='CSS-ML')
coeftest(model_1216_ml_css)
residual.analysis(model = model_1216_ml_css)

# AIC and BIC
sc.AIC=AIC(model_718_ml, model_818_ml, model_418_ml, model_419_ml, model_518_ml, model_12112_ml, 
           model_12110_ml, model_1216_ml)
sc.BIC=AIC(model_718_ml, model_818_ml, model_418_ml, model_419_ml, model_518_ml, model_12112_ml, 
           model_12110_ml, model_1216_ml,  k = log(length(diff.employment)))
sc.BIC
sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "aic") 

# choose best model: using ARIMA(12, 1, 10)

m519Residuals = model_12110_ml$residuals
abs.res = abs(m519Residuals)
sq.res = m519Residuals^2


par(mfrow=c(1,2))
acf(abs.res, main="The ACF plot for absolute residual series")
pacf(abs.res, main="The PACF plot for absolute residual series")
par(mfrow=c(1,1))

# -------SEASONAL ARIMA-------------------
# Helper function ---------------------------------------------------------------------

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

#Work with data from 2000 to 2020 only
#so models would run faster and 20 years are still
#good enough for forecasting
subset = window(employment_TS,2000,2020)
subset

sqrt_data <- sqrt(subset)

sqrt_data
plot(sqrt_data,ylab='Time series plot of BC transformed 
     yearly average employment numbers.')

# acf and pacf plots
par(mfrow=c(1,2))
show_ACF_PACF(sqrt_data,maxlag=48,name="Squared-rooted subset series")
par(mfrow=c(1,1))


# DIFFERENCING

diff_sqrt = diff(sqrt_data)
par(mfrow=c(1,1))
plot(diff_sqrt,type='o',ylab='Average employment numbers', main = "Time series plot of the first differenced
     yearly average employment numbers.")
abline(h=0)
adf.test(diff_sqrt) # stationary
pp.test(diff_sqrt) # stationary
kpss.test(diff_sqrt) # stationary
#All tests indicate the series at first difference is stationary

#SEASONAL ARIMA
#Start with a baseline model with D=1 and d=1
m1 = Arima(subset,order=c(0,0,0),seasonal=list(order=c(0,1,0), period=12))
res.m1 = residuals(m1);  
par(mfrow=c(1,1))
plot(res.m1,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
show_ACF_PACF(res.m1, maxlag = 48, name="Residuals of model 1")
#From the ACF and PACF => P=2 and Q=

#
m2 = Arima(subset,order=c(0,0,0),seasonal=list(order=c(2,1,3), period=12),method="CSS")
res.m2 = residuals(m2)
par(mfrow=c(1,1))
plot(res.m2,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
show_ACF_PACF(res.m2, maxlag = 48,name="Residuals of model 2")
adf.test(res.m2) #Residuals are stationary

#Use sqrt_data for less variance
#Add differencing
m3 = Arima(sqrt_data,order=c(0,1,0),seasonal=list(order=c(2,1,3), period=12))
res.m3 = residuals(m3);  
par(mfrow=c(1,1))
plot(res.m3,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
show_ACF_PACF(res.m3, maxlag = 48,name="Residuals of model 3")
#No trend or seasonality left, move forward
adf.test(res.m3)


#From ACF and PACF of m3
#p = 3 and q = 5
#Determined by significant bars between 0 and 1
m4 = Arima(sqrt_data,order=c(3,1,5),
           seasonal=list(order=c(2,1,3), period=12),
           method="CSS")
res.m4 = residuals(m4);  
par(mfrow=c(1,1))
plot(res.m4,xlab='Time',ylab='Residuals',main="Time series plot of the residuals")
show_ACF_PACF(res.m4, maxlag = 48,name="Residuals of model 4")
# SARIMA(3,1,5)x(2,1,3)_12

eacf(res.m3) 
# The tentative models are specified as 
# SARIMA(0,1,1)x(2,1,3)_12
# SARIMA(0,1,2)x(2,1,3)_12
# SARIMA(1,1,1)x(2,1,3)_12
# SARIMA(1,1,2)x(2,1,3)_12

par(mfrow=c(1,1))
bic_table = armasubsets(y=res.m3,nar=10,nma=10,y.name='p',ar.method='ols')
plot(bic_table)
# SARIMA(8,1,0)x(2,1,3)_12
# SARIMA(1,1,0)x(2,1,3)_12

#Final sets of possible models
# SARIMA(3,1,5)x(2,1,3)_12
# SARIMA(0,1,1)x(2,1,3)_12
# SARIMA(0,1,2)x(2,1,3)_12
# SARIMA(1,1,1)x(2,1,3)_12
# SARIMA(1,1,2)x(2,1,3)_12
# SARIMA(8,1,0)x(2,1,3)_12
# SARIMA(1,1,0)x(2,1,3)_12

#Define a function to fit models quickly
fit_SARIMA <- function(series,orders) {
  model_css = Arima(series,order=orders,
                    seasonal=list(order=c(2,1,3),period=12),
                    method='CSS')
  coef_css = coeftest(model_css)
  print("CSS Method")
  print(coef_css)
  
  model_ml = Arima(series,order=orders,
                   seasonal=list(order=c(2,1,3),period=12),
                   method='ML')
  coef_ml = coeftest(model_ml)
  print("ML Method")
  print(coef_ml)
  return(list(css=model_css,ml=model_ml))
}


# SARIMA(0,1,1)x(2,1,3)_12
m011 = fit_SARIMA(sqrt_data,orders=c(0,1,1))
#ML - all variables are significant
residual.analysis(model = m011$css) #Did not capture the info well

# SARIMA(0,1,2)x(5,1,1)_12
m012=fit_SARIMA(sqrt_data,orders=c(0,1,2))
residual.analysis(model = m012$css)  #Did not capture the info well

# SARIMA(1,1,1)x(5,1,1)_12
m111 = fit_SARIMA(sqrt_data,orders=c(1,1,1))
#ML - NaNs in sma2 and sma3
residual.analysis(model = m111$css) #Did not capture the info well

#CSSML
m111mlcss = Arima(sqrt_data,order=c(1,1,1),
                  seasonal=list(order=c(2,1,3), period=12),method = "CSS-ML")
coeftest(m111mlcss) 
residual.analysis(model = m111mlcss)

# SARIMA(1,1,2)x(5,1,1)_12
m112 = fit_SARIMA(sqrt_data,orders=c(1,1,2))
residual.analysis(model = m112$css) #Did not capture the info well
m112mlcss = Arima(sqrt_data,order=c(1,1,2),
                  seasonal=list(order=c(2,1,3), period=12),method = "CSS-ML")
coeftest(m112mlcss)
residual.analysis(model = m112mlcss) #Did not capture the info well

x11()
# SARIMA(3,1,5)x(2,1,3)_12
m315 = fit_SARIMA(sqrt_data,orders=c(3,1,5))
residual.analysis(model = m315$ml) #pretty residuals

# SARIMA(8,1,0)x(5,1,1)_12
m810 = fit_SARIMA(sqrt_data,orders=c(8,1,0))
#Not all variables are significant
residual.analysis(model = m810$css) #pretty good residuals

# SARIMA(1,1,0)x(2,1,3)_12
m110 =fit_SARIMA(sqrt_data,orders=c(1,1,0)) #ML - all vars are significant
residual.analysis(model = m110$css) # Did not capture info well


sc.AIC = AIC(m011$ml,m012$ml,m111$ml,m112$ml,m810$ml,m315$ml,m110$ml)
sc.BIC = BIC(m011$ml,m012$ml,m111$ml,m112$ml,m810$ml,m315$ml,m110$ml)
sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "bic")

#Residuals of best models from BIC still have
#significant autocorrelation

#Best models from AIC captured the information well
#and nothing valuable is in their residuals

#Try overfitted model
m415 = fit_SARIMA(sqrt_data,orders=c(4,1,5)) #Most vars are significant in CSS
residual.analysis(m415$ml) #Pretty good residuals


m316 = fit_SARIMA(sqrt_data,orders=c(3,1,6)) #About half of the vars are significant
residual.analysis(m315$ml) #Pretty good residuals

#ERROR METRICS
m011css <- accuracy(m011$css)[1:7]
m012css <- accuracy(m012$css)[1:7]
m111css <- accuracy(m111$css)[1:7]
m112css <- accuracy(m112$css)[1:7]
m315css <- accuracy(m315$css)[1:7]
m810css <- accuracy(m810$css)[1:7]
m110css <- accuracy(m110$css)[1:7]
m415css <- accuracy(m415$css)[1:7]
m315css <- accuracy(m315$css)[1:7]
df.Smodels <- data.frame(
  rbind(m011css, m012css, m111css, m112css, m315css, m810css, m110css,
        m415css,m315css)
)
colnames(df.Smodels) <- c("ME", "RMSE", "MAE", "MPE", "MAPE", 
                          "MASE", "ACF1")
rownames(df.Smodels) <- c("SARIMA(0,1,1)x(2,1,3)_12", "SARIMA(0,1,2)x(2,1,3)_12", "SARIMA(1,1,1)x(2,1,3)_12", 
                          "SARIMA(1,1,2)x(2,1,3)_12", "SARIMA(3,1,5)x(2,1,3)_12", "SARIMA(8,1,0)x(2,1,3)_12",
                          "SARIMA(1,1,0)x(2,1,3)_12",  "SARIMA(4,1,5)x(2,1,3)_12", "SARIMA(3,1,6)x(2,1,3)_12")
round(df.Smodels,  digits = 3)


sc.AIC = AIC(m011$ml,m012$ml,m111$ml,m112$ml,
             m810$ml,m315$ml,m110$ml,m415$ml,m316$ml)
sc.BIC = BIC(m011$ml,m012$ml,m111$ml,m112$ml,
             m810$ml,m315$ml,m110$ml,m415$ml,m316$ml)
sort.score(sc.AIC, score = "aic")
sort.score(sc.BIC, score = "bic")

#AIC picks m316 over 315

#FORECASTING
forecastML = forecast(m316$ml,  h = 10)
forecastML
plot(forecastML)

forecastCSS = forecast(m316$css,  h = 10)
forecastCSS
plot(forecastCSS) #forecast my CSS model has a sligher smaller interval
