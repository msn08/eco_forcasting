##### ********************************* ####
##### Part 0: Packages, functions-files #### 
##### ********************************* ####
#*Packackes to install:
#*    To read the data
#*install.packages("readr")

#*    For the plots
#*install.packages("dplyr")
#*install.packages("ggplot2")
#*install.packages("scales")

#*    For the timeseries
#*install.packages("xts")
#*install.packages("tseries") # adf.test, na.omit() 
#*install.packages("forecast") # autoplot(),snaive(), auto.arima(),Arima(),
#* # forecast(), tsCV()
#*install.packages("MLmetrics") # MAPE()
#*install.packages("astsa") # acf2()


#*   Functions from "functions_project6"
#*   Packages needed:
#*install.packages("sandwich") # vcovHAC
#*install.packages("lmtest") # coeftest
#*install.packages("car") # linearHypothesis
#*install.packages("data.table") # data.table()

#*   Functions:
#* pasta_dataframe_exo() 
#* pasta_plot 
#* pasta_forecast_matrix()
#* matrix_errors_exo()
#* mz_test()
#* tau_test()
#* ic_renne()
#* ic_saron()
#* get_forecast_evaluation_exo() and all its dependecies
#* auto_forecast() and all it's dependencies

##### ********************************* ####
##### Part 1: Exploratory Data Analysis ####
##### ********************************* ####
library(readr)
data <- read_csv("https://raw.githubusercontent.com/msn08/eco_forcasting/main/data2/COVID19Cases_geoRegion.csv")

#****  NATIONAL LEVEL
#****
library(dplyr)
library(ggplot2)
library(scales)
by_ch <- subset(data, geoRegion == "CH")

plot_ch <- ggplot(by_ch,aes(x=as.Date(datum), y = entries,fill=entries))+
  scale_fill_distiller(palette="Reds",direction=1)+
  geom_bar(width=0.7,stat="identity")+
  theme_minimal()

plot_ch + labs(title="Total number of COVID-19 cases in Switzerland")+
  xlab("Days")+
  ylab("Total Cases per day")


#****  TIMESERIES 
#****

#****  SUBSETTING TO MATCH VACC DATA 
#****
by_ch$datum <- as.Date(by_ch$datum)

#same dates 
by_ch_date <- subset(by_ch, datum >= as.Date("2020-12-21") & 
                       datum <= as.Date("2021-11-24"))


library(xts)
covid_ts <- xts(by_ch_date$entries,order.by=as.Date(by_ch_date$datum, "%Y%m/%d/"),frequency=7)

#Plot Data and Deterministic Trend
plot(covid_ts,type="l",main="Total COVID-19 Cases per Day",ylab="Cases",xlab="Data")

library(tseries)
covid_ts = ts(by_ch_date$entries , frequency = 7, start = c(2020,355))
plot(covid_ts,type="l",main="Total COVID-19 Cases per Day",ylab="Cases",xlab="Data")
abline(reg=lm(covid_ts~time(covid_ts)), col="red")

#We see we from 4000 cases to 6, so we'll remove the last one
to_remove <- length(covid_ts)
covid_ts[to_remove] <- NA
covid_ts <- na.remove(covid_ts)
plot(covid_ts,type="l",main="Total COVID-19 Cases per Day",ylab="Cases",xlab="Data")
abline(reg=lm(covid_ts~time(covid_ts)), col="red")

#We take the trend-cycle component
compo.covid <- decompose(covid_ts)
plot(compo.covid)

covid_ts <- compo.covid$trend
#Getting Vaccine Data
vacc <- read.csv("https://raw.githubusercontent.com/msn08/eco_forcasting/main/data2/COVID19VaccDosesAdministered.csv")
vacc_by_ch <- subset(vacc, geoRegion == "CH")

#****  NATIONAL LEVEL
#****
vacc_by_ch <- subset(vacc, geoRegion == "CH")

plot_ch <- ggplot(vacc_by_ch,aes(x=as.Date(date), y = entries,fill=entries))+
  scale_fill_distiller(palette="Reds",direction=1)+
  geom_bar(width=0.7,stat="identity")+
  theme_minimal()

plot_ch + labs(title="Total number of COVID-19 Vaccinations in Switzerland")+
  xlab("Days")+
  ylab("Total Vaccinations per day")

vacc_ts = ts(vacc_by_ch$entries , frequency = 7, start = c(2020,355))
to_remove <- length(vacc_ts)
vacc_ts[to_remove] <- NA
vacc_ts <- na.remove(vacc_ts)

#Both in one graph
library(forecast)
case_vacc <- cbind(covid_ts,vacc_ts)
colnames(case_vacc) <- c("Daily COVID Cases","Daily Vaccination Doses")
autoplot(case_vacc,facets=TRUE)+
  xlab("Days")+ ylab("")+
  ggtitle("Daily COVID cases and Vaccinaition Doses")


##### ********************************* ####
##### Part 2: ARIMAX Model and Analysis ####
##### ********************************* ####
#**** Check for stationarity
#**** Original Data
library(tseries)
covid_ts <- na.omit(covid_ts)
adf.test(covid_ts, alternative = "stationary")
# not stationary

#**** d=1
covid_ts_diff<- diff(covid_ts)
covid_ts_diff <- na.omit(covid_ts_diff)
# not stationary, alternative = "stationary")
# It is stationnary

#**** d=2
covid_ts_diff<- diff(diff(covid_ts))
covid_ts_diff <- na.omit(covid_ts_diff)
adf.test(covid_ts_diff, alternative = "stationary")
# It is stationnary

#**** d=3
covid_ts_diff<- diff(diff(diff(covid_ts)))
covid_ts_diff <- na.omit(covid_ts_diff)
adf.test(covid_ts_diff, alternative = "stationary")
# It is stationnary


#**** Modelling
library(forecast)
library(MLmetrics)
##### d= 0, Baseline to compare later

covid_ts_diff<- covid_ts
#**** Splitting into training data, and test data
#* This allows us to do a real out-of-example evaluation
case_vacc <- cbind(covid_ts_diff,vacc_ts)
colnames(case_vacc) <- c("Daily COVID Cases","Daily Vaccination Doses")
total_length <- length(covid_ts)
periods <- 30
end_train <- total_length - periods
start_test <- end_train + 1

#### Final two datasets
in_sample  <- window(case_vacc[0:end_train,1])
out_sample <- window(case_vacc[start_test:total_length,1])


#### Model and Forecasting
# The lags for the exogenous variable. Previous analysis shows that
# Adlag = 3 is best. This is consistent with the medical data 
vacc_lagged <- cbind(
  AdLag0 = case_vacc[,"Daily Vaccination Doses"],
  AdLag1 = stats::lag(case_vacc[,"Daily Vaccination Doses"],-7),
  AdLag2 = stats::lag(case_vacc[,"Daily Vaccination Doses"],-14),
  AdLag3 = stats::lag(case_vacc[,"Daily Vaccination Doses"],-21),
  AdLag4 = stats::lag(case_vacc[,"Daily Vaccination Doses"],-28),
  AdLag5 = stats::lag(case_vacc[,"Daily Vaccination Doses"],-35),
  AdLag6 = stats::lag(case_vacc[,"Daily Vaccination Doses"],-42),
  AdLag7 = stats::lag(case_vacc[,"Daily Vaccination Doses"],-49),
  AdLag8 = stats::lag(case_vacc[,"Daily Vaccination Doses"],-56),
  AdLag9 = stats::lag(case_vacc[,"Daily Vaccination Doses"],-63))%>%
  head(NROW(case_vacc))

#*****   Random walk, no exogenous variable, d=0
naive = snaive(in_sample, h=periods)
plot(naive)
summary(naive)
qqnorm(naive$residuals)
MAPE(naive$mean, out_sample) * 100
#out of sample MAPE: 44.06573
#Ilias's computer:
#out of sample MAPE: 

#*****   AR1, exogenous variable, d=0
ar1 <- c(1,0,0)
for(i in 1:10){
  print("LAG")
  print(i)
  table <- get_forecast_evaluation_exo(in_sample,vacc_lagged[0:end_train,1:i],ar1,periods,
                                       out_sample,vacc_lagged[start_test:total_length,1:i])
  print(table)
}
#out of sample MAPE: 42.33782
# i=1 best, meaning one lag (Seven days after vaccination)
#Ilias's computer:
#out of sample MAPE: 

#*****   ARIMAX MODEL SELECTION
#* Need to change d in auto_forecast() manually

for(i in 1:10){
  print("LAG")
  print(i)
  table <- auto_forecast(in_sample,vacc_lagged[0:end_train,1:i],periods,
                         out_sample,vacc_lagged[start_test:total_length,1:i])
}
#**** d= 3 yields best results, this was found on a previous analysis
# Model selected ARIMA (4,3,2)

#*****  ARIMA (4,3,2)
# This model was selected by auto_forecast()( by lowest MAPE)  when d = 3
fit_basic1 <- Arima(in_sample, order=c(4,3,2),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE 13.08691
#Ilias's computer:
#out of sample MAPE: 


#*****  ARIMA (5,3,0)
# This model was selected by auto.arima()  when d = 3
fit_basic1 <- Arima(in_sample, order=c(5,3,0),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE: 37.03596
#Ilias's computer:
#out of sample MAPE: 

#*****  ARIMA (12,3,0)
# This model was selected by ic_saron()( by lowest AIC)  when d = 3
fit_basic1 <- Arima(in_sample, order=c(12,3,0),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE: 25.0754
#Ilias's computer:
#out of sample MAPE: 

#*****  ARIMA (6,3,6)
# This model was selected by ic_renne()( by lowest AIC)  when d = 3
fit_basic1 <- Arima(in_sample, order=c(6,3,6),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE: 15.63536
#Ilias's computer:
#out of sample MAPE: 

# We see that our model has the best real out of sample performance

#### Other stationary models
# d = 2 models

#*****  ARIMA (4,2,0)
# This model was selected by auto.arima()  when d = 2
fit_basic1 <- Arima(in_sample, order=c(4,2,0),xreg=vacc_lagged[0:end_train,3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE: 15.04313
#Ilias's computer:
#out of sample MAPE: 

#*****  ARIMA (12,2,0)
# This model was selected by ic_saron()( by lowest AIC)  when d = 2
fit_basic1 <- Arima(in_sample, order=c(12,2,0),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE: 28.66112
#Ilias's computer:
#out of sample MAPE: 

#*****  ARIMA (6,2,6)
# This model was selected by ic_renne()( by lowest AIC)  when d = 2
fit_basic1 <- Arima(in_sample, order=c(6,2,6),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE: 19.87583
#Ilias's computer:
#out of sample MAPE: 

#We see that our model still has the best real out of sample performance

###### PSEUDO OUT OF SAMPLE EVALUATION
#*****  ARIMA (4,3,2)
covid_ts_diff<- diff(diff(diff(covid_ts)))

#Create initial pasta_dataframe
pasta_covid<- data.frame(date=1:length(covid_ts_diff), u=covid_ts_diff)
params <-c(4,0,2)
length(covid_ts_diff)

xreg_in <- vacc_lagged[,1:3]
length(xreg_in)
#Create pasta_dataframe with forecasts
new_pasta_covid_exo <- pasta_dataframe_exo(pasta_covid,330,4,30,params,xreg_in)

pasta_plot(new_pasta_covid_exo,330,4,30)

#Get the forecasts into a matrix
new_matrix_for_exo <- pasta_forecast_matrix(new_pasta_covid_exo,330,4,30)

#Matrix of forecast errors
params <- c(4,0,2)
y <- new_pasta_covid_exo[1:338,2]#to  match the length
exoregressor <- vacc_lagged[,3]
length(y)
length(exoregressor)

# Matrix of forecast errors, but model needs to be speficied by hand!!!!
# If another model is test, go to functions_project6.r to change the parameters
error_new <- matrix_errors_exo(y,30,exoregressor)
#error_new <- matrix_errors(y,params,30)

#Evaluation
## Mincer-Zarnowitz (MZ) Regression
mz_test(new_pasta_covid_exo[,2],new_matrix_for_exo,330,4,30)
#bad 

## Tau Test
tau_test(error_new,30)
#[1,]  -24.986127 -2.4797956 0.01367335 #biased
#[7,]  -23.288015 -1.7960556 0.07347203
#[14,]  -85.080195 -1.2553836 0.21031974
#[30,] -121.233515 -1.4571110 0.14619442
# unbiased

#Ilias's computer:

## test on 1 step ahead forecast errors follows a MA(1)

acf(error_new[,1], na.action = na.pass) 
pacf(error_new[,1], na.action = na.pass)
# Not optimal


acf(error_new[,7], na.action = na.pass) 
pacf(error_new[,7], na.action = na.pass)
# Not optimal

acf(error_new[,15], na.action = na.pass) 
pacf(error_new[,15], na.action = na.pass)
# optimal

acf(error_new[,30], na.action = na.pass) 
pacf(error_new[,30], na.action = na.pass)
# optimal

#Our forcast for 30 days is unbiased and weakly optimal 

##### ********************************* ####
##### Part 3: Decisions based on models ####
##### ********************************* ####
#We will make evaluate the following policies:
# Policy 1: No change is made, vaccinations continues it's course
# Policy 2: 2G policies which increase vaccination by 20%
# Policy 3: Full lockdown for unvaccinated, which increases vaccination by 40%
# Policy 4: Lift all restrictions which reduces vaccinations by 50%

covid_ts_diff<- diff(diff(diff(covid_ts)))
y <- covid_ts_diff[1:338]#to  match the length
exoregressor <- vacc_lagged[,3]
length(y)


### Base line, one week before forecast
mean_baseline <- mean(in_sample[297:303])


#*****  Policy 1
lenght_exo <- length(exoregressor)
start_exo <- lenght_exo - 30
policy_1_vax <- exoregressor[start_exo:lenght_exo]
fit_basic1 <- Arima(covid_ts, order=c(6,3,6),xreg=exoregressor)
pred_0<-forecast(fit_basic1, h=30,xreg = policy_1_vax )$mean
plot(pred_0)
policy1_full <-pred_1

#*****  Policy 2
policy_2_vax <- policy_1_vax * 2
fit_basic1 <- Arima(covid_ts, order=c(6,3,6),xreg=exoregressor)
pred_1<-forecast(fit_basic1, h=30,xreg = policy_2_vax )$mean
plot(pred_1)
policy1_full <-pred_1



policy_1_vax <- exoregressor
fit_basic1 <- arima(covid_ts, order=c(12,2,0),xreg=exoregressor)
pred_2<-predict(fit_basic1, n.ahead=30,newxreg = policy_1_vax )$pred
plot(pred_2)
policy2_full <-pred_2

policy_2_vax <- policy_1_vax * 2
fit_basic1 <- arima(covid_ts, order=c(12,2,0),xreg=exoregressor)
pred_3<-predict(fit_basic1, n.ahead=30,newxreg = policy_2_vax )$pred
plot(pred_3)
policy2_full <-pred_2
