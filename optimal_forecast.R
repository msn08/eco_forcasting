####### *********  ********************  ********* #######
####### *********  ECONOMIC FORECASTING  ********* ####### 
####### *********  ********************  ********* #######

#####Project by Ilias and Marco Antonio

#####Libraries
#install.packages("readr")
#install.packages("dplyr")
#install.packages("tseries") 
#install.packages("astsa")
#install.packages("forecast")
#install.packages("MLmetrics")

library(readr)
data <- read_csv("https://raw.githubusercontent.com/msn08/eco_forcasting/main/data2/COVID19Cases_geoRegion.csv")
#Description: Daily record timelines by geoRegion for cases.

####### *********  EXPLORATORY DATA ANALYSIS  ********* ####### 

colnames(data)
#descriptions of columns: https://www.covid19.admin.ch/api/data/documentation/models/sources-definitions-dailyincomingdata.md
#*Columns we care about:
#**** getRegion:  
#*Cantons and top level Units CHFL, CH or FL for aggregated records.
#*Those need to be removed

#****  entries: 
#*Absolute number of occurrences for this day

#**** mean7d: 
#*Rolling 7 day average (-3 days/ +3 days from the current day) of the "entries" values.
#*mean7d could be used later, but not there
library(dplyr)
library(ggplot2)
library(scales)

#****  CANTON LEVEL
#****
by_region_day <- subset(data, 
                        geoRegion != "CH" & geoRegion != "FL" & geoRegion != "CHFL")

total_canton <- aggregate(by_region_day$entries, 
                          by=list(Category=by_region_day$geoRegion), FUN=sum)


plot_canton <- ggplot(total_canton,aes(x=Category, y = x,fill=x))+
  scale_fill_distiller(name = "Total\nCases",palette="Reds",
                       direction=1,breaks = c(50000, 149000))+
  geom_bar(width=0.7, stat = "identity",color="grey")+
  scale_y_continuous(labels = unit_format(unit = "K", scale = 1e-3))

plot_canton + labs(title="Total number of COVID-19 cases by Canton")+
  xlab("Cantons")+
  ylab("Total Cases in thousands")

#****  NATIONAL LEVEL
#****
by_ch <- subset(data, geoRegion == "CH")

plot_ch <- ggplot(by_ch,aes(x=datum, y = entries,fill=entries))+
  scale_fill_distiller(palette="Reds",direction=1)+
  geom_bar(width=0.7,stat="identity")+
  theme_minimal()

plot_ch + labs(title="Total number of COVID-19 cases in Switzerland")+
  xlab("Days")+
  ylab("Total Cases per day")

#****  TIMESERIES 
#****
library(xts)
covid_ts <- xts(by_ch$entries,order.by=as.Date(by_ch$datum, "%Y%m/%d/"),frequency=7)

#Plot Data and Deterministic Trend
plot(covid_ts,type="l",main="Total CODVID-19 Cases per Day",ylab="Cases",xlab="Data")

covid_ts = ts(by_ch$entries , frequency = 365, start = c(2020,55))
plot(covid_ts,type="l",main="Total CODVID-19 Cases per Day",ylab="Cases",xlab="Data")
abline(reg=lm(covid_ts~time(covid_ts)), col="red")


####### *********  TIMESERIES ANALYSIS  ********* ####### 

### STEP 0: GETTING STATIONARY DATA

#****  Vizualize Trends and Seasonality
#****
covid_ts = ts(by_ch$entries , frequency = 7, start = c(2020,55))
compo.covid <- decompose(covid_ts)
plot(compo.covid)

#****  Select trend component for our analysis
#****
covid_trend <- compo.covid$trend
compo.covid_trend <- decompose(covid_trend)
plot(compo.covid_trend)

#****  Remove Seasonality
#****
covid_trend_se <- covid_trend-compo.covid_trend$seasonal
plot(covid_trend_se, type = "l")
abline(reg=lm(covid_trend_se~time(covid_trend_se)), col="blue")
#Still upward trend -> not stationary 

#****  Detrend by difference
#****
covid_trend_diff <- diff(covid_trend_se)
plot(covid_trend_diff, type = "l", main = "Detrended and deseasonalized COVID-19 CASES", ylab = "Daily Cases", xlab="Time, in days")
abline(reg=lm(covid_trend_diff~time(covid_trend_diff)), col="blue")

#**** Check for stationarity
#****
library(tseries)
covid_trend_diff <- na.omit(covid_trend_diff)
adf.test(covid_trend_diff, alternative = "stationary")

### STEP 1: SPLITTING THE DATA
###
library(MLmetrics)

#**** Create samples
#****
total_length = length(covid_trend_diff)
periods = 30
end_train = total_length - periods
start_test = end_train + 1
in_sample=window(covid_trend_diff[0:end_train])
out_sample=window(covid_trend_diff[start_test:total_length])


### STEP 2: Model Selection
###
library(forecast)

#**** ACF and PACF Analysis
#****
Pacf(in_sample, lag.max = 90)
Acf(in_sample, lag.max = 90)

library(astsa)
acf2(in_sample, main="ACF and PACF")

#**** BIC and AIC
#****
##### d = 1
ic_renne(diff(in_sample))
#AIC: p = 4, q = 4
#BIC: p = 4, q = 2


ic_saron(diff(in_sample))
#AIC: p = 12, q = 2
#BIC: p = 9, q = 1

##### d = 0
ic_renne(in_sample)
#AIC: p = 5, q = 5
#BIC: p = 5, q = 5

ic_saron(in_sample)
#AIC: p = 11, q = 2
#BIC: p = 9, q = 2

#we see significant changes given d, specially with ic_renne
#maybe investigate more

#**** Automatic Model Selection
#****
model_auto <- auto.arima(in_sample,method="ML")
summary(model_auto)
#p = 3, q = 1

model_auto_diff <- auto.arima(diff(in_sample),method="ML")
summary(model_auto_diff)
#p = 2, q = 0

### STEP 2.1: Model Selection, it's complicated
### We have lots of disagreement regarding the models
### We will fit all of them and then compared them 
### Papers also recommend to do +/- 1 in the neighbourg of auto.arima

auto_table <- auto_neighbour(in_sample,3,0,1)
#We get new model from bic: arima(3,0,0)


auto_table_diff <- auto_neighbour(in_sample,2,1,0)
#We get the same model!
#First time aic and bic agree, those pieces of shit
#arima(2,1,0)

#* d=0
#renne_m1 = arima(5,0,5)
#saron_m1 = arima(11,0,2)
#saron_m2 = arima(9,0,2)
#model_auto = arima(3,0,1)
#model_auto_table = arima(3,0,0)
renne_m1 = c(5,0,5)
saron_m1 = c(11,0,2)
saron_m2 = c(9,0,2)
model_auto = c(3,0,1)
model_auto_table = c(3,0,0)
ar1 = c(1,0,0)


#C'EST MOCHE, FAUT AUTOMATISER FLEMMARD
table_compare <- get_ic(in_sample,5,0,5)
table_compare2 <- get_ic(in_sample,11,0,2)
total_table <- cbind(table_compare,table_compare2)
new_table <- get_ic(in_sample,9,0,2)
total_table <- cbind(total_table,new_table)
new_table <- get_ic(in_sample,3,0,1)
total_table <- cbind(total_table,new_table)
new_table <- get_ic(in_sample,3,0,0)
total_table <- cbind(total_table,new_table)
get_neighbours(total_table)

#We get two models
#From AIC: saron_m1 = arima(11,0,2)
#Fromb BIC: saron_m2 = arima(9,0,2)

#* d=1


#renne_diff_m1 = arima(4,1,4)
#renne_diff_m2 = arima(4,1,2)
#saron_diff_m1 = arima(12,1,2)
#saron_diff_m2 = arima(9,1,1)
#model_auto_diff = arima(2,1,0)

renne_diff_m1 = c(4,1,4)
renne_diff_m2 = c(4,1,2)
saron_diff_m1 = c(12,1,2)
saron_diff_m2 = c(9,1,1)
model_auto_diff = c(2,1,0)
ar1 = c(1,0,0)

table_compare <- get_ic(in_sample,4,1,4)
table_compare2 <- get_ic(in_sample,4,1,2)
total_table <- cbind(table_compare,table_compare2)
new_table <- get_ic(in_sample,12,1,2)
total_table <- cbind(total_table,new_table)
new_table <- get_ic(in_sample,9,1,1)
total_table <- cbind(total_table,new_table)
new_table <- get_ic(in_sample,2,1,1)
total_table <- cbind(total_table,new_table)
get_neighbours(total_table)

#We get two models
#From AIC: renne_diff_m2 = arima(11,0,2)
#Fromb BIC: saron_diff_m2 = arima(9,0,2)


### STEP 3: Forecasting Selection
#Random walk
naive = snaive(in_sample, h=periods)
plot(naive)
summary(naive)
qqnorm(naive$residuals)
MAPE(naive$mean, out_sample) * 100

#* d = 0
list_models = list(renne_m1,saron_m1,saron_m2,
                   model_auto,model_auto_table,ar1)
for (model in list_models){
  get_forecast_evaluation(in_sample,model,periods,out_sample)
}

#* d = 1
list_models = list(renne_diff_m1,renne_diff_m2,saron_diff_m1,
                   saron_diff_m2,model_auto_diff,ar1)
for (model in list_models){
  get_forecast_evaluation(in_sample,model,periods,out_sample)
}

#We get the model that beats both an AR(1) and a RW
# arima(2,1,0)
#technically (2,2,0) since we differentiated once already
fit_basic1 <- Arima(in_sample, order=c(1,0,0))
forecast_1<-forecast(fit_basic1, h=periods)
plot(forecast_1)
plot(forecast_1$residuals)
acf(forecast_1$residuals)
pacf(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
MAPE(forecast_1$mean, out_sample) * 100

