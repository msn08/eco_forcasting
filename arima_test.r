####### *********  ********************  ********* #######
####### *********  ECONOMIC FORECASTING  ********* ####### 
####### *********  ********************  ********* #######

#####Project by Ilias and Marco Antonio

#####Libraries
#install.packages("readr")
library(readr)
#install.packages("dplyr")
library(dplyr)
library(ggplot2)
library(scales)
library(xts)
#install.packages("tseries") 
#library(tseries)
#tseries ca marche pas chez moi faut que je vois avec mon autre ordi
#install.packages("astsa")
library(astsa)
#install.packages("forecast")
library(forecast)
#forecast ca marche pas chez moi faut que je vois avec mon autre ordi

####### *********  ARIMA BUT NO EXO YET  ********* ####### 

#DATA EXPLORATION
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


by_region_day <- subset(data, 
                        geoRegion != "CH" & geoRegion != "FL" & geoRegion != "CHFL")
View(total_canton)
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
covid_ts_se <- covid_ts - compo.covid$seasonal
plot(covid_ts_se, type = "l")
abline(reg=lm(covid_ts_se~time(covid_ts_se)), col="blue")

#Detrend
y <- covid_ts_se
trend <- 1:length(y)
eq <- lm(y~trend)
y <- eq$residuals
plot(y, type="l")
covid_ts_se <- y

plot(covid_ts_se, type = "l", main = "Detrended and deseasonalized COVID-19 CASES", ylab = "Daily Cases", xlab="Time, in days")
abline(reg=lm(covid_ts_se~time(covid_ts_se)), col="blue")

#Check for stationairy
#ici on a besoin de tseries
adf.test(covid_ts_se, alternative = "stationary")

#ACF and PACF
acf2(covid_ts_se, main="ACF and PACF")


#model estimation
#ici on a besoin de forecast
auto.arima(covid_ts_se, seasonal = FALSE, trace = TRUE) 

#BIC,AIC,HQ
y <- covid_ts_se
plot(y,type="l")

T <- length(y)
max.p <- 5
max.q <- 5
matrix.AIC <- matrix(NaN,max.p+1,max.q+1)
matrix.BIC <- matrix(NaN,max.p+1,max.q+1)
matrix.HQ  <- matrix(NaN,max.p+1,max.q+1)
for(p in 0:max.p){
  for(q in 0:max.q){
    fit.arima <- arima(y,order=c(p,1,q),method="ML")
    matrix.AIC[p+1,q+1] <- -2*fit.arima$loglik/T + (p+q+2)*2/T
    matrix.BIC[p+1,q+1] <- -2*fit.arima$loglik/T + (p+q+2)*log(T)/T
    matrix.HQ[p+1,q+1]  <- -2*fit.arima$loglik/T + (p+q+2)*2*log(log(T))/T
  }
}
print(min(matrix.AIC))
print(matrix.AIC)
#ARiMA(6,1,6)
print(min(matrix.BIC))
print(matrix.BIC)
#ARiMA(6,1,6)
print(min(matrix.HQ))
print(matrix.HQ)
#ARiMA(6,1,6)

#*PROBLEME:
#*LES TROIS DONNE LES VALEURS MAX POUR LES PARAMETRES
#*LE ACF ET PACF N'AIDENT PAS
#*L'IDEAL SERAIT AUTO.ARIMA MAIS CA MARCHE PAS SUR MON ORDI MAINTENANT. 
#*NORMALEMENT LE MOINS DE PARAMETRE POSSIBLE, LE MIEUX, JE PENSE Y A TROP
