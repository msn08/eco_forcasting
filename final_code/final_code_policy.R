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
#*install.packages("lubridate") #to get dates in autoplot correct

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
#to_remove <- length(covid_ts)
#covid_ts[to_remove] <- NA
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

library(tseries)
vacc_ts = ts(vacc_by_ch$entries , frequency = 7, start = c(2020,355))
compo.vacc <- decompose(vacc_ts)
plot(compo.vacc)
vacc_ts <- compo.vacc$trend
#to_remove <- length(vacc_ts)
#vacc_ts[to_remove] <- NA
vacc_ts2 = ts(vacc_ts , frequency = 365, start = c(2020,355))
covid_ts2 = ts(covid_ts , frequency = 365, start = c(2020,355))
#Both in one graph

library(forecast)
library(lubridate)
case_vacc2 <- cbind(covid_ts2,vacc_ts2)
colnames(case_vacc2) <- c("Daily COVID Cases","Daily Vaccination Doses")
my_date_transform <- function(x) {format(date_decimal(x), "%m/%d/%y")}

autoplot(case_vacc2, facets=TRUE)+
  xlab("Days")+ ylab("")+
  ggtitle("Daily COVID cases and Vaccination Doses")+
  scale_x_continuous(labels = my_date_transform)


#Both in one graph, wrong dates, but you have to run or the rest won't run
library(forecast)
case_vacc <- cbind(covid_ts,vacc_ts)
colnames(case_vacc) <- c("Daily COVID Cases","Daily Vaccination Doses")
autoplot(case_vacc,facets=TRUE)+
  xlab("Days")+ ylab("")+
  ggtitle("Daily COVID cases and Vaccination Doses")

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
#out of sample MAPE: 46.13684

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
#out of sample MAPE: 39.08318
# i=10

#Ilias's computer:
#out of sample MAPE: 


#*****   ARIMAX MODEL SELECTION (IN SAMPLE PERFORMANCE)
#auto.arima() from forecast library
auto.arima(in_sample,xreg=vacc_lagged[0:end_train,1:i])
#Selection: ARIMA(2,2,0)

#ic_renne from Pr. Renne's course "Macroeconomics"
#d=3
diff_ts <- diff(diff(diff(in_sample)))
ic_renne(diff_ts)
#AIC Selection: ARIMA (4,3,4)
#BIC Selection: ARIMA(2,3,0)

#d=2
diff_ts <- diff(diff(in_sample))
ic_renne(diff_ts)
#AIC Selection: ARIMA (3,2,4)
#BIC Selection: ARIMA(1,2,3)

#ic_saron from Pr. Grobéty's course "Economic Forecasting for Decision Making"
#d=3
diff_ts <- diff(diff(diff(in_sample)))
ic_saron(diff_ts)
#AIC Selection: ARIMA (8,3,2)
#BIC Selection: ARIMA(0,3,2)

#d=2
diff_ts <- diff(diff(in_sample))
ic_saron(diff_ts)
#AIC Selection: ARIMA (8,2,1)
#BIC Selection: ARIMA(3,2,1)

#*****   ARIMAX MODEL SELECTION (OUT OF SAMPLE PERFORMANCE)
#* Need to change d in auto_forecast() manually

for(i in 1:10){
  print("LAG")
  print(i)
  table <- auto_forecast(in_sample,vacc_lagged[0:end_train,1:i],periods,
                         out_sample,vacc_lagged[start_test:total_length,1:i])
}

#**** d= 3 yields best results, this was found on a previous analysis
# Model selected ARIMA:
#1_3_1 i=3
#6_3_1 i=3
#*****  ARIMA (1,3,1)
fit_basic1 <- Arima(in_sample, order=c(1,3,1),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE 9.20857
#Ilias's computer:
#out of sample MAPE: 

#*****  ARIMA (6,3,1)
fit_basic1 <- Arima(in_sample, order=c(6,3,1),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE 8.877753
#Ilias's computer: 

#auto_arima
fit_basic1 <- Arima(in_sample, order=c(2,2,0),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE 13.21535

#ic_renne

fit_basic1 <- Arima(in_sample, order=c(4,3,4),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE 10.44777


fit_basic1 <- Arima(in_sample, order=c(2,3,0),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE 12.43695

fit_basic1 <- Arima(in_sample, order=c(3,2,4),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE: 20.40329


fit_basic1 <- Arima(in_sample, order=c(1,2,3),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE: 12.66341

#ic_saron 
#aic
fit_basic1 <- Arima(in_sample, order=c(8,3,2),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE: 10.98751
#bic
fit_basic1 <- Arima(in_sample, order=c(0,3,2),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE: 9.648655

#aic
fit_basic1 <- Arima(in_sample, order=c(8,2,1),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE: 15.06924

#bic
fit_basic1 <- Arima(in_sample, order=c(3,2,1),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
plot(forecast_1)
plot(forecast_1$residuals)
library(astsa)
acf2(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
print("MAPE:")
print(MAPE(forecast_1$mean, out_sample) * 100)
#out of sample MAPE:13.03318

### OUT OF ALL THE MODELS, ARIMA(6,3,1) is best 

######
###### PSEUDO OUT OF SAMPLE EVALUATION
######
#*****  ARIMA (1,3,1)
covid_ts_diff<- diff(diff(diff(in_sample)))

#Create initial pasta_dataframe
pasta_covid<- data.frame(date=1:length(covid_ts_diff), u=covid_ts_diff)
params <-c(6,0,1)
length(covid_ts_diff)

xreg_in <- vacc_lagged[,1:3]
length(xreg_in)
#Create pasta_dataframe with forecasts
new_pasta_covid_exo <- pasta_dataframe_exo(pasta_covid,300,7,30,params,xreg_in)

pasta_plot(new_pasta_covid_exo,300,7,30)

#Get the forecasts into a matrix
new_matrix_for_exo <- pasta_forecast_matrix(new_pasta_covid_exo,300,7,30)

#Matrix of forecast errors
params <- c(6,0,1)
y <- new_pasta_covid_exo[1:339,2]#to  match the length
exoregressor <- vacc_lagged[,3]
length(y)
length(exoregressor)

# Matrix of forecast errors, but model needs to be speficied by hand!!!!
# If another model is test, go to functions_project6.r to change the parameters
error_new <- matrix_errors_exo(y,30,exoregressor)
#error_new <- matrix_errors(y,params,30)

#Evaluation
## Mincer-Zarnowitz (MZ) Regression
mz_test(new_pasta_covid_exo[,2],new_matrix_for_exo,300,7,30)
#bad 

## Tau Test
tau_test(error_new,30)
#  [1,]  -4.6119928 -1.6859671 0.09292007
#  [7,]  -1.9982827 -0.9077041 0.36483512
# [14,]  -3.2704107 -1.0176302 0.30977840
# [28,] -12.0185085 -1.1805535 0.23889767
# unbiased

#Ilias's computer:

## test on 1 step ahead forecast errors follows a MA(1)

acf(error_new[,1], na.action = na.pass) 
pacf(error_new[,1], na.action = na.pass)
# Not optimal


acf(error_new[,7], na.action = na.pass) 
pacf(error_new[,7], na.action = na.pass)
#  optimal

acf(error_new[,14], na.action = na.pass) 
pacf(error_new[,14], na.action = na.pass)
# optimal

acf(error_new[,28], na.action = na.pass) 
pacf(error_new[,28], na.action = na.pass)
# optimal

#DM Test

library(forecast)
#*****  ARIMA (6,3,1)
fit_basic1 <- Arima(in_sample, order=c(6,3,1),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
model <- forecast_1$mean

#*****  ARIMA (1,0,0)
fit_basic1 <- Arima(in_sample, order=c(1,0,0),xreg=vacc_lagged[0:end_train,1:3])
forecast_1<-forecast(fit_basic1, h=periods,xreg= vacc_lagged[start_test:total_length,1:3])
ar1 <- forecast_1$mean

dm.test(model,ar1,h=28)

#Our forcast for 30 days is unbiased and weakly optimal 

##### ********************************* ####
##### Part 3: Decisions based on models ####
##### ********************************* ####

#We will make evaluate the following policies:
# Policy 1: No change is made, vaccinations continues it's course
# Policy 2: A policy which increases vaccination by 38%
# Policy 3: Full lockdown for unvaccinated, which double vaccination 
# Policy 4: Lift all restrictions which reduces vaccinations by 50%
# Policy 5: Triples the current rate, in line with data

exoregressor <- vacc_lagged[1:333,3]

length(exoregressor)


### Base line, one week before forecast
swiss_pop <- 8670300
vax_base <- (sum(na.omit(exoregressor))/2)
vax_perc <- vax_base/swiss_pop


#*****  Policy 1

exoregressor <- vacc_lagged[0:end_train,3]

policy_1_vax <- vacc_lagged[start_test:total_length,3]
fit_basic1 <- Arima(in_sample, order=c(6,3,1),xreg=exoregressor)
pred_1<-forecast(fit_basic1, h=60,xreg = policy_1_vax )$mean
pred_1
policy1_full <-pred_1

avg_in_sample <- mean(in_sample[273:303])
avg_pol1 <- mean(policy1_full)
(avg_pol1/avg_in_sample)*100#how much it goes up

sum_pol1 <- sum(policy_1_vax)
policy1_week1 <- mean(policy1_full[1:7])
policy1_week2 <- mean(policy1_full[7:14])
policy1_week3 <- mean(policy1_full[15:21])
policy1_week4 <- mean(policy1_full[22:30])

#*****  Policy 2
policy_2_vax <- policy_1_vax * 1.38
fit_basic1 <- Arima(in_sample, order=c(6,3,1),xreg=exoregressor)
pred_2<-forecast(fit_basic1, h=60,xreg = policy_2_vax )$mean
pred_2
policy2_full <-pred_2
avg_pol2 <- mean(policy2_full)
(avg_pol2/avg_in_sample)*100
policy2_week1 <- mean(policy2_full[1:7])
policy2_week2 <- mean(policy2_full[7:14])
policy2_week3 <- mean(policy2_full[15:21])
policy2_week4 <- mean(policy2_full[22:30])
sum_pol2 <- sum(policy_2_vax)
pol2_perc <- (sum_pol2 + vax_base)/swiss_pop

#*****  Policy 3
policy_3_vax <- policy_1_vax * 2
fit_basic1 <- Arima(in_sample, order=c(6,3,1),xreg=exoregressor)
pred_3<-forecast(fit_basic1, h=60,xreg = policy_3_vax )$mean
pred_3
policy3_full <-pred_3
avg_pol3 <- mean(policy3_full)
(avg_pol3/avg_in_sample)*100
policy3_week1 <- mean(policy3_full[1:7])
policy3_week2 <- mean(policy3_full[7:14])
policy3_week3 <- mean(policy3_full[15:21])
policy3_week4 <- mean(policy3_full[22:30])
sum_pol3 <- sum(policy_3_vax)
pol3_perc <- (sum_pol3 + vax_base)/swiss_pop

#*****  Policy 4
policy_4_vax <- policy_1_vax * 0.5
fit_basic1 <- Arima(in_sample, order=c(6,3,1),xreg=exoregressor)
pred_4<-forecast(fit_basic1, h=60,xreg = policy_4_vax )$mean
pred_4
policy4_full <-pred_4
avg_pol4 <- mean(policy4_full)
(avg_pol4/avg_in_sample)*100
policy4_week1 <- mean(policy4_full[1:7])
policy4_week2 <- mean(policy4_full[7:14])
policy4_week3 <- mean(policy4_full[15:21])
policy4_week4 <- mean(policy4_full[22:30])
sum_pol4 <- sum(policy_4_vax)
pol4_perc <- (sum_pol4 + vax_base)/swiss_pop

#*****  Policy 5
policy_5_vax <- policy_1_vax * 3.31
fit_basic1 <- Arima(in_sample, order=c(6,3,1),xreg=exoregressor)
pred_5<-forecast(fit_basic1, h=60,xreg = policy_5_vax )$mean
pred_5
policy5_full <-pred_5
avg_pol5 <- mean(policy5_full)
(avg_pol5/avg_in_sample)*100
sum(policy1_full)-sum(policy5_full)
policy5_week1 <- mean(policy5_full[1:7])
policy5_week2 <- mean(policy5_full[7:14])
policy5_week3 <- mean(policy5_full[15:21])
policy5_week4 <- mean(policy5_full[22:30])
sum_pol5 <- sum(policy_5_vax)
pol5_perc <- (sum_pol5 + vax_base)/swiss_pop

### Graphs
pol2 <- c((1-(policy2_week1/policy1_week1))*100,(1-(policy2_week2/policy1_week2))*100,
          (1-(policy2_week3/policy1_week3))*100,(1-(policy2_week4/policy1_week4))*100)
pol3 <- c((1-(policy3_week1/policy1_week1))*100,(1-(policy3_week2/policy1_week2))*100,
          (1-(policy3_week3/policy1_week3))*100,(1-(policy3_week4/policy1_week4))*100)
pol4 <- c((1-(policy4_week1/policy1_week1))*100,(1-(policy4_week2/policy1_week2))*100,
          (1-(policy4_week3/policy1_week3))*100,(1-(policy4_week4/policy1_week4))*100)
pol5 <- c((1-(policy5_week1/policy1_week1))*100,(1-(policy5_week2/policy1_week2))*100,
          (1-(policy5_week3/policy1_week3))*100,(1-(policy5_week4/policy1_week4))*100)


plot(pol2,ylim=c(-5,20),type="l",col="black",xlab="Weeks",
     ylab="Reduction of Daily Cases in Percent",main="4 Policies Compared")
lines(pol3,col="green")
lines(pol4,col="red")
lines(pol5,col="blue")

pop_perc <- c(vax_perc,pol2_perc,pol3_perc,pol4_perc,pol5_perc)
labels <- c("Baseline","Policy-2","Policy-3","Policy-4","Policy-5")

df <- data.frame(Policies=labels,Coverage=pop_perc)
p<-ggplot(data=df, aes(x=Policies, y=Coverage)) +
  geom_bar(stat="identity")
p

plot(in_sample,ylim=c(900,3000),xlim=c(283,330),type="l",col="black",xlab="Days",
     ylab="Daily Cases",main="Future Situation")
lines(policy1_full,col="red")
lines(policy5_full,col="blue")


#THIS IS HOW I KNOW IT IS CORRECT

###### Literature review on the impact of vaccines

#The effect of mandatory COVID-19 certificates on vaccine uptake: synthetic-control modelling of six countries
#Prof Melinda C Mills, PhD, Tobias Rüttenauer

# Mandatory COVID-19 certification could increase vaccine uptake, but interpretation and transferability of findings need to be considered in the context of pre-existing levels of vaccine uptake and hesitancy, eligibility changes, and the pandemic trajectory.

#In Switzerland, approximately 1 month before the introduction of certificates, vaccination levels were lower than in control countries and daily doses were above average briefly before the intervention (appendix p 26). Again, vaccination rates exceeded the average of the control countries for 40 days after introduction of certification (47 380 [95% CI 9870–78 627] doses per million population; appendix p 26). Results based on different estimation methods support the long-lasting upward shift after certification (appendix p 23).

#Switzerland first introduced some access restrictions (events >1000 participants and nightclubs), and later extended restrictions to general situations and activities (appendix p 7). These earlier access restrictions appear to have only affected vaccine uptake in those younger than 20 years, while uptake did not change among older groups (figure 4). Extending these restrictions to more general activities continued to affect uptake in those younger than 20 years, but also uptake among older age groups (30–39 years and 40–49 years; figure 4).

#Appendix: table a.5.1: anticipation 153'152 more vax, 20 days before certificates

##Notes
#This is about 3.31 the compare to the one month prior. (153152/mean(vacc_ts[106:126]))
#Why compared it to one month prior?
#Because of the anticipations, we see from articles at the time that this is indeed the case
#https://www.rts.ch/info/suisse/12181531-le-certificat-covid-sera-realise-par-ladministration-federale-dici-fin-juin.html
#https://www.rts.ch/info/suisse/12167210-pharmasuisse-et-la-fmh-lancent-leur-propre-certificat-covid19.html
#These are from 20 days before the announcement. The article says 20 days before the implementation of the policy
#Since we don't the know the actual effects, for simplicity we take these numbers.


#The Effect of Vaccination Rates on the Infection of COVID-19 under the Vaccination Rate below the Herd Immunity Threshold
#Yi-Tui Chen

#However, the infection rate after vaccination showed two trends. One is an inverted U-shaped trend, and the other is an L-shaped trend. For those countries with an inverted U-shaped trend, the infection rate begins to decline when the vaccination rate reaches 1.46-50.91 doses per 100 people.

#This means that higher vaccination rates in the UK and the USA led to lower infection rates during the whole analysis period, and the theoretical peak takes places when vaccination programs start.

#Only when the accumulated vaccination rate reaches a certain level does partial protection of herd immunity more or less take place, and the infection of disease is reduced

#It suggests that without using a mask, a 50% effective vaccine will not suppress infection at low vaccination coverage rate, but an 80% effective vaccine requires a 48–78% vaccination coverage rate, and a 100% effective vaccine requires 33–58% vaccination coverage to curb the spread of the COVID-19 pandemic. In the case of a mask usage rate of 50%, a 50% effective vaccine requires a 55–94% vaccination coverage rate, while an 80% effective vaccine only requires a 32–57% vaccination coverage rate, and a 100% effective vaccine requires a 24–46% vaccination coverage rate to suppress the spread of the pandemic

#This paper suggests that vaccines are not a panacea in solving the pandemic of COVID-19 if the new variant provides a significant impact on transmissibility, severity, and/or immunity. If individuals currently vaccinated are proven to lack resistance to new variants, new vaccines should be developed and, thus, the pandemic is unlikely to be stopped in a short time.

##Notes
#The peaked was reached in june
#https://www.covid19.admin.ch/en/vaccination/doses
# required at least 78% in switzerland, we are at 68.26%

#Which policy will be able to reach which level?

#Policy 2: "Vaccine week policy"
#-> drive by to get vaccines
#increases about 38%
#https://www.swissinfo.ch/eng/swiss--vaccine-week--considered-only-partial-success/47113188

