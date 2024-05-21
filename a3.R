rm(list=ls())
library(TSA)
library(tseries)
library(forecast)
library(lmtest)
library(dplyr)

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

# Calculate the Z-scores
z_scores <- scale(employment_TS)

# Identify outliers (e.g., Z-score > 3 or < -3)
outliers_z <- which(abs(z_scores) > 3)
outliers_z

# acf and pacf plots
par(mfrow=c(1,2))
acf(employment_TS,lag.max = 98, main ="ACF plot for the Global Temperature Anomalies series.")
pacf(employment_TS, main ="PACF plot for the Global Temperature Anomalies series.")
par(mfrow=c(1,1))


# qq plot 
qqnorm(employment_TS)
qqline(employment_TS, col = 2)
shapiro.test(employment_TS)


# adf test (weak indicator of non-stationary)
adf.test(employment_TS, alternative = c("stationary")) # Use default value of k

# pp test
pp.test(employment_TS)
kpss.test(employment_TS)

summary(employment_TS)


BC <- BoxCox.ar(employment_TS, lambda = seq(0.5, 1, 0.01))
BC$ci
lambda <- BC$lambda[which(max(BC$loglike) == BC$loglike)]
lambda # lambda = 0.5 corresponds to square root transformation?
# BC.employment = (employment_TS^lambda-1)/lambda

# plot(BC.employment,type='o',ylab='Time series plot of BC transformed 
# yearly average unemployment numbers.')

sqrt_data <- sqrt(employment_TS)

plot(sqrt_data,type='o',ylab='Time series plot of BC transformed 
     yearly average unemployment numbers.')
# A small positive constant
log_diff_employment <- log(sqrt_data + 2.7)
min_value <- min(sqrt_data)
min_value

# qq plot
qqnorm(log_diff_employment)
qqline(log_diff_employment, col = 2)
shapiro.test(log_diff_employment)

par(mfrow=c(1,1))
qqnorm(sqrt_data)
qqline(sqrt_data, col = 2)
shapiro.test(sqrt_data)

par(mfrow=c(1,2))
acf(sqrt_data, main='ACF plot of employment series.')
pacf(sqrt_data, main='PACF plot of employment series.')
par(mfrow=c(1,1))

diff.employment = diff(sqrt_data)

par(mfrow=c(1,1))
plot(diff.employment,type='o',ylab='Average employment numbers', main = "Time series plot of the first differenced
     yearly average employment numbers.")
# There is still some trend in there series
adf.test(diff.employment) # stationary
pp.test(diff.employment) # stationary
kpss.test(diff.employment) # non stationary

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


par(mfrow=c(1,2))
acf(diff.employment, main='ACF plot of unemployment series.')
pacf(diff.employment, main='PACF plot of unemployment series.')
par(mfrow=c(1,1))


# qq plot 
qqnorm(diff.employment)
qqline(diff.employment, col = 2)
shapiro.test(diff.employment)

eacf(diff.employment)
res = armasubsets(y=diff.employment , nar=5 , nma=5, y.name='p', ar.method='ols')
plot(res)



mean_diff_employment <- mean(diff.employment, na.rm = TRUE)
sd_diff_employment <- sd(diff.employment, na.rm = TRUE)

# Standardize the data
standardized_diff_employment <- (sqrt_data - mean_diff_employment) / sd_diff_employment

# Plot the standardized data
plot(sqrt_data, type = 'o', ylab = 'Standardized Differenced Employment', main = "Time series plot of standardized first differenced yearly average employment numbers")

# Check normality of the standardized data
par(mfrow=c(1,1))
qqnorm(standardized_diff_employment)
qqline(standardized_diff_employment, col = 2)
shapiro.test(standardized_diff_employment)
# # check if there are non positive values
# min_value <- min(diff.employment)
# min_value
# 
#  # A small positive constant
# log_diff_employment <- log(diff.employment + 2.7)
# min_value <- min(diff.employment)
# min_value
# 
# # qq plot 
# qqnorm(log_diff_employment)
# qqline(log_diff_employment, col = 2)
# shapiro.test(log_diff_employment)



# # Normalize function
# normalize <- function(x) {
#   return((x - min(x)) / (max(x) - min(x)))
# }
# 
# # Normalize the time series data
# normalized_ts <- normalize(diff.employment)
# 
# diff.employment
# # qq plot 
# qqnorm(normalized_ts)
# qqline(normalized_ts, col = 2)
# shapiro.test(normalized_ts)




# too many significant lags, difference again?
# 
# diff.employment2 = diff(diff(sqrt_data))
# 
# par(mfrow=c(1,1))
# plot(diff.employment2,type='o',ylab='Average employment numbers', main = "Time series plot of the second differenced
#      yearly average employment numbers.")
# # There is still some trend in there series
# adf.test(diff.employment2) # stationary
# pp.test(diff.employment2) # stationary
# kpss.test(diff.employment2) # non stationary
# 
# 
# par(mfrow=c(1,2))
# acf(diff.employment, main='ACF plot of unemployment series.')
# pacf(diff.employment, main='PACF plot of unemployment series.')
# par(mfrow=c(1,1))
