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


####### *********  ARIMA BUT NO EXO YET  ********* ####### 

#DATA EXPLORATION
library(readr)
data <- read_csv("https://raw.githubusercontent.com/msn08/eco_forcasting/main/data/COVID19Cases_geoRegion.csv")
#Description: Daily record timelines by geoRegion for cases.

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

##** PREMIERE INTEGORRATION:
##* VU LA DIFFERENCE ENTRE LES CANTONS, EST-CE QU'ON DEVRAIT FAIRE DES REGRESSION
##* PAR CANTON, OU SUR TOUTE LA SUISSE?
##* LA JE PARS SUR TOUTE LA SUISSE

by_ch <- subset(data, geoRegion == "CH")

plot_ch <- ggplot(by_ch,aes(x=datum, y = entries,fill=entries))+
  scale_fill_distiller(palette="Reds",direction=1)+
  geom_bar(width=0.7,stat="identity")+
  theme_minimal()

plot_ch + labs(title="Total number of COVID-19 cases in Switzerland")+
  xlab("Days")+
  ylab("Total Cases per day")

#TIME SERIES
#get start time
library(xts)
covid_ts <- xts(by_ch$entries,order.by=as.Date(by_ch$datum, "%Y%m/%d/"),frequency=7)

#Plot Data and Deterministic Trend
plot(covid_ts,type="l",main="Total CODVID-19 Cases per Day",ylab="Cases",xlab="Datw")
#abline doesn't work with xts
covid_ts = ts(by_ch$entries , frequency = 365, start = c(2020,55))
plot(covid_ts,type="l",main="Total CODVID-19 Cases per Day",ylab="Cases",xlab="Datw")
abline(reg=lm(covid_ts~time(covid_ts)), col="red")

#Vizualize Trends and Seasonality

#NOTE: THE DECOMPOSE FUNCTION NEEDS TO HAVE AT LEAST 2 "CYCLE"
#SINCE OUR DATA IS DAILY, THE FREQUENCY SHOULD BE 365
#BUT, THAT DOESN'T COMPLETE THE 2 CYCLES AND THUS DECOMPOSE DOENS'T WORK
#BY HAVING FREQUENCY = 7, IT COUNT'S UP TO 2120, SO INSTEAD OF A YEAR 
#IT COUNTS 100 YEARS, AND THE FUNCTION DECOMPOSE WORKS.
#THEN WE JUST NEED TO READJUST THE DATA POINTS
covid_ts = ts(by_ch$entries , frequency = 7, start = c(2020,55))
compo.covid <- decompose(covid_ts)
plot(compo.covid)


#Remove Sesonality
covid_trend <- compo.covid$trend
compo.covid_trend <- decompose(covid_trend)
plot(compo.covid_trend)

covid_trend_se <- covid_trend-compo.covid_trend$seasonal
plot(covid_trend_se, type = "l")
abline(reg=lm(covid_trend_se~time(covid_trend_se)), col="blue")


#Detrend by Difference
covid_trend_diff <- diff(covid_trend_se)

plot(covid_trend_diff, type = "l", main = "Detrended and deseasonalized COVID-19 CASES", ylab = "Daily Cases", xlab="Time, in days")
abline(reg=lm(covid_trend_diff~time(covid_trend_diff)), col="blue")

#Check for stationairy
library(tseries)
covid_trend_diff <- na.omit(covid_trend_diff)
adf.test(covid_trend_diff, alternative = "stationary")

#ACF and PACF
library(astsa)
acf2(covid_trend_diff, main="ACF and PACF")


#MODEL SELECTION
library(forecast)
auto.arima(covid_trend_diff, seasonal = FALSE, ic="bic",trace = TRUE, method ="ML") 
# Best model: ARIMA(3,0,0)           with zero mean     

#BIC,AIC,HQ

ic_renne(covid_trend_diff) #(2,0,3)
ic_saron(covid_trend_diff) #(9,0)
BIC(covid_trend_diff)
### We have three candidates
### let's use sarima() for the selction 
#From auto arima
sarima(covid_trend_diff, 3,0,0)
#BIC: 10.78311
#Visual inspection:
#acf of residuals: bad
#awful q-q plot
#Ljung box: most values under the line
#This one is trash

#From ic_renne
sarima(covid_trend_diff, 2,0,3)
#BIC:10.77488
#Visual inspection:
#acf of residuals: bad but better
#awful q-q plot
#Ljung box: most values under the line
#This one is also trash

#From ic_saron
sarima(covid_trend_diff, 9,0,0)
#BIC: 10.7729
#Visual inspection:
#acf of residuals: much better, only two values cross the line
#awful q-q plot but better
#Ljung box: most values above the line
#This one is not realy good but the best as of now




####### *********  FORECASTING  ********* ####### 

library(MLmetrics)


#Create samples
in_sample=window(covid_trend_diff[0:528])
out_sample=window(covid_trend_diff[529:619])

periods = length(out_sample)

#Random walk
naive = snaive(in_sample, h=periods)
plot(naive)
summary(naive)
qqnorm(naive$residuals)
MAPE(naive$mean, out_sample) * 100
#

#AR1 forecast
fit_basic1 <- arima(in_sample, order=c(1,0,0))
forecast_1<-forecast(fit_basic1, h=periods)
plot(forecast_1)
plot(forecast_1$residuals)
acf(forecast_1$residuals)
pacf(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
MAPE(forecast_1$mean, out_sample) * 100


#Auto.arima
plot(in_sample, type = "l", main = "Detrended and deseasonalized COVID-19 CASES", ylab = "Daily Cases", xlab="Time, in days")
fit_basic1<- auto.arima(in_sample)
#here auto.arima chooses ARIMA(3,0,1)
#Forecast to see preliminary results 
forecast_1<-forecast(fit_basic1, h=periods)
plot(forecast_1)
plot(forecast_1$residuals)
acf(forecast_1$residuals)
pacf(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
MAPE(forecast_1$mean, out_sample) * 100

#We do not beat AR1 but we beat RW


#ic_saron
fit_basic1 <- arima(in_sample, order=c(9,0,0))
forecast_1<-forecast(fit_basic1, h=periods)
plot(forecast_1)
plot(forecast_1$residuals)
acf(forecast_1$residuals)
pacf(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
MAPE(forecast_1$mean, out_sample) * 100

#We do not beat AR1 but we beat RW

#ic_renne
fit_basic1 <- arima(in_sample, order=c(2,0,3))
forecast_1<-forecast(fit_basic1, h=periods)
plot(forecast_1)
plot(forecast_1$residuals)
acf(forecast_1$residuals)
pacf(forecast_1$residuals)
summary(fit_basic1)
qqnorm(forecast_1$residuals)
MAPE(forecast_1$mean, out_sample) * 100

#We do not beat AR1 but we beat RW

