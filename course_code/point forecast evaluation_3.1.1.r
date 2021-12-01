rm(list=ls())

library(RCurl)#for downloading the data
library(readxl)
library(foreign)
library(lmtest)
library(sandwich)
library(car)
library(data.table)







employment_rate <- read_excel("employment_china.xlsx", col_types = c("date", rep("numeric", 3)))
employment_rate<-employment_rate[-c(1:5),]
employment_rate<-employment_rate[-c(1:43),]
employment_rate<-employment_rate[1:74,]
names(employment_rate)[4]<-"u"
setDT(employment_rate)

#AR(3) model

AR3<-arima(employment_rate$u, order = c(3,0,0))

#Spaghetti dataframe

Spaghetti_dataframe_u<-data.frame(date=employment_rate$...1, u=employment_rate$u)

arima(Spaghetti_dataframe_u[1:4,2], order = c(3,0,0),method="ML")

for (c in 1:71) {

  model.fit <- arima(Spaghetti_dataframe_u[1:(3+c), 2], order = c(3, 0, 0), optim.method="Nelder-Mead" )
  Spaghetti_dataframe_u[(4+c):(3+c+12),2+c] <- predict(model.fit, n.ahead = 12)$pred
  
}

#Spaghetti graph
plot(Spaghetti_dataframe_u[1:(74+12), 2], type="l", xaxt="n", xlab="", ylab="unemployment" )
for (i in 1:71){
  
  lines(Spaghetti_dataframe_u[(1):(74+12),(2+i)], type = "l", col="lightblue")
  
}
axis(side=1, at=c(6,26,46,66,86), labels = c(2004, 2009, 2014, 2019, 2024) )

# Transformation of the Spaghetti dataframe into matrix

m_forecast<-matrix(nrow = 71, ncol = 12)

m_forecast<-matrix_dataframe(m_forecast, Spaghetti_dataframe_u, 71, 3, 12)

# Matrix of forecast errors

m_error<-matrix(nrow = 71, ncol = 12)

for (c in 1:71){
  
  m_error[c,]<-Spaghetti_dataframe_u[(4+c):(15+c),2]-m_forecast[c,]
  
}

# unbiasedness and efficiency

## Mincer-Zarnowitz (MZ)

modelMZ<-c()
vector_result_MZ<-c()

x<-lm(Spaghetti_dataframe_u$u[(4+1):(62+1)]~m_forecast[1:59,1]) 
linearHypothesis(x,c("(Intercept)=0", "m_forecast[1:59, 1]=1"))$`Pr(>F)`[2]

for (i in 1:12) {
  
  modelMZ[[i]]<- lm(Spaghetti_dataframe_u$u[(4+i):(62+i)]~m_forecast[1:59,i])
  vector_result_MZ[i]<-linearHypothesis(modelMZ[[i]], c("(Intercept)=0", "m_forecast[1:59, i]=1"))$`Pr(>F)`[2]
  
}

##test: H_0: tao=0
m_tao_result<-matrix(nrow = 12,ncol = 3)
colnames(m_tao_result)<-c("coefficient","t-value", "p-value")

for (i in 1:12){

  model<-lm(m_error[,i]~ 1)  
  result <-coeftest(model, vcov = vcovHAC(model))
  m_tao_result[i,1]<-result[1]
  m_tao_result[i,2]<-result[3]
  m_tao_result[i,3]<-result[4]

}

##test on available information

modelINF<-c()
vector_result_INF<-c()

for (i in 1:12) {
  modelINF[[i]]<-lm(m_error[1:59,i]~ Spaghetti_dataframe_u$u[(1+i):(59+i)]+Spaghetti_dataframe_u$u[(2+i):(60+i)]+Spaghetti_dataframe_u$u[(3+i):(61+i)] )
  vector_result_INF[i]<-linearHypothesis(modelINF[[i]], c("Spaghetti_dataframe_u$u[(1 + i):(59 + i)]=0", "Spaghetti_dataframe_u$u[(2 + i):(60 + i)]=0", "Spaghetti_dataframe_u$u[(3 + i):(61 + i)]=0"))$`Pr(>F)`[2]
}

## test on 1 step ahead forecast errors follows a MA(1)

acf(m_error[,1], na.action = na.pass) #MA(1)
pacf(m_error[,1], na.action = na.pass)

## test on h steps ahead forecast errors follows at most MA(h-1)

for (i in 2:12) {
  
  acf(m_error[,i], na.action = na.pass) 
  pacf(m_error[,i], na.action = na.pass)
  
}
