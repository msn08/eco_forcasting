####Functions for AIC & BIC

#From Pr Renne
ic_renne <- function(data_ts){
  y <- data_ts
  plot(y,type="l")
  
  T <- length(y)
  max.p <- 5
  max.q <- 5
  matrix.AIC <- matrix(NaN,max.p+1,max.q+1)
  matrix.BIC <- matrix(NaN,max.p+1,max.q+1)
  matrix.HQ  <- matrix(NaN,max.p+1,max.q+1)
  for(p in 0:max.p){
    for(q in 0:max.q){
      fit.arima <- arima(y,order=c(p,0,q),method="ML")
      matrix.AIC[p+1,q+1] <- -2*fit.arima$loglik/T + (p+q+2)*2/T
      matrix.BIC[p+1,q+1] <- -2*fit.arima$loglik/T + (p+q+2)*log(T)/T
      matrix.HQ[p+1,q+1]  <- -2*fit.arima$loglik/T + (p+q+2)*2*log(log(T))/T
    }
  }
  print("MIN AIC")
  print(min(matrix.AIC))
  print(matrix.AIC)
  #ARiMA(5,0,5)
  print("MIN BIC")
  print(min(matrix.BIC))
  print(matrix.BIC)
  #ARiMA(2,0,3)
  print("MIN HQ")
  print(min(matrix.HQ))
  print(matrix.HQ)
  #ARiMA(5,0,5)
}

#From 2.2.1_GM-5.14.1_saron_Switzerland.r
library(data.table)
ic_saron <- function(data_column){
  # ARIMA selection
  ## loop computes the AIC and BIC of each ARIMA with a maximum of AR lag = 12
  ## and a maximum of MA lag = 2
  ic.saron <- list('AIC' = data.table(), 'BIC' = data.table())
  for (ar.lag in 0:12) {
    arma.stat <- rep(0, 6)
    for (ma.lag in 0:2) {
      arma.fit <- arima(data_column, order = c(ar.lag, 0, ma.lag),method="ML")
      arma.stat[ma.lag + 1] <- arma.fit$aic
      arma.stat[ma.lag + 4] <- -2 * arma.fit$loglik +
        (ar.lag + ma.lag) * log(214)
    }
    ic.saron$AIC <- rbindlist(list(ic.saron$AIC, data.table(t(arma.stat[1:3]))))
    ic.saron$BIC <- rbindlist(list(ic.saron$BIC, data.table(t(arma.stat[4:6]))))
  }
  setnames(ic.saron$AIC, c('MA0', 'MA1', 'MA2'))
  ic.saron$AIC[, AR := (0:12)]
  setnames(ic.saron$BIC, c('MA0', 'MA1', 'MA2'))
  ic.saron$BIC[, AR := (0:12)]
  
  arima.fit <- list()
  print("MIN AIC")
  print(min(ic.saron$AIC[,1:3]))
  print("MIN BIC")
  print(min(ic.saron$BIC[,1:3]))
  print(ic.saron)
}

