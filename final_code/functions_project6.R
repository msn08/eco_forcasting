library(forecast)
library(MLmetrics)#MAPE
library(car)#linearHypothesis
library(lmtest)#coeftest
library(sandwich)#vcovHAC
######## GET INFO CRITERA

#compute model
get_ic <- function(data,list_parameters){
  model = Arima(data, order = c(list_parameters))		
  model_aic <- model$aicc
  model_bic <- model$bic
  names_ic <- c("aic","bic")
  list_ic <- c(model_aic,model_bic)
  name_col <- paste(c(list_parameters), collapse = "_")
  final_list <- as.data.frame(list_ic,row.names=names_ic)
  colnames(final_list) =name_col
  return(final_list)
}

#makes table with min IC
get_table_ic <- function(data_table){
  print(data_table)
  print("Min AIC:")
  print(min(data_table[1,]))
  print("Min BIC:")
  print(min(data_table[2,]))
}

#does +/- 1 for p and q, and applied get_ic() and get_table_ic()
auto_neighbour <- function(data, list_parameters){
  p <- list_parameters[1]
  d <- list_parameters[2]
  q <- list_parameters[3]
  if(p == 0){
    i <- 0
  }
  else{
    i <- -1
  }
  table <- get_ic(data,list_parameters)
  while(i <= 1){
    new_p <- p + i
    if(q == 0){
      j <- 0
    }
    else{
      j <- -1
    }
    while(j <= 1){
      new_q <- q + j
      new_list_parameters <- c(new_p,d,new_q)
      new_table <- get_ic(data,new_list_parameters)
      table <- cbind(table,new_table)
      j <- j + 1
    }
    i <- i + 1
  }

  get_table_ic(table)
  return(table)
}

#similar to auto_neighbor but with pre determined models
#always starts with ar1
auto_compare <- function(data,list_models){
  ar1 <- c(1,0,0)
  new_table <- get_ic(data,ar1)
  for (i in 1:length(list_models)){
    model <- unlist(list_models[i])
    table <- get_ic(data,model)
    new_table <- cbind(new_table,table)
  }
  get_table_ic(new_table)
  return(new_table)
}


######## GET INFO CRITERA + MAPE FORCASTING PART


### NO EXO
get_forecast_evaluation <- function(data,list_parameters,periods_forecast,out_of_sample){
  #get model and forecast
  model = Arima(data, order = c(list_parameters))	
  forecast_1<-forecast(model, h=periods_forecast)
  #ic
  model_aic <- model$aicc
  model_bic <- model$bic
  #mape
  model_mape <- MAPE(forecast_1$mean, out_sample) * 100
  #make table
  names_ic <- c("aic","bic","mape")
  list_ic <- c(model_aic,model_bic,model_mape)
  name_col <- paste(c(list_parameters), collapse = "_")
  final_list <- as.data.frame(list_ic,row.names=names_ic)
  colnames(final_list) =name_col
  return(final_list)
}

get_table_min <- function(data_table){
  print(data_table)
  print("Min AIC:")
  print(min(data_table[1,]))
  print("Min BIC:")
  print(min(data_table[2,]))
  print("Min MAPE:")
  print(min(data_table[3,]))
}

auto_compare_mape <- function(data,list_of_models,periods_forecast,out_of_sample){
  ar1 <- c(1,0,0)
  new_table <- get_forecast_evaluation(data,ar1,periods_forecast,out_of_sample)
  for (i in 1:length(list_of_models)){
    model <- unlist(list_of_models[i])
    table <- get_forecast_evaluation(data,model,periods_forecast,out_of_sample)
    new_table <- cbind(new_table,table)
  }
  get_table_min(new_table)
  return(new_table)
}

## EXO 
get_forecast_evaluation_exo <- function (data,xreg_in,list_parameters,periods_forecast,
                                         out_of_sample,xreg_out) {
  model <- Arima(data, order=c(list_parameters),xreg=xreg_in)
  forecast_1<-forecast(model, h=periods_forecast,xreg= xreg_out)
  #ic
  model_aic <- model$aicc
  model_bic <- model$bic
  #mape
  model_mape <- MAPE(forecast_1$mean, out_sample) * 100
  #make table
  names_ic <- c("aic","bic","mape")
  list_ic <- c(model_aic,model_bic,model_mape)
  name_col <- paste(c(list_parameters), collapse = "_")
  final_list <- as.data.frame(list_ic,row.names=names_ic)
  colnames(final_list) =name_col
  return(final_list)
}

auto_compare_exo_mape <- function (data,xreg_in,list_of_models,periods_forecast,
                                    out_of_sample,xreg_out){
  ar1 <- c(1,0,0)
  new_table <- get_forecast_evaluation_exo(data,xreg_in,ar1,
                                           periods_forecast,
                                           out_of_sample,xreg_out)
  for (i in 1:length(list_of_models)){
    model <- unlist(list_of_models[i])
    table <-  get_forecast_evaluation_exo(data,xreg_in,model,
                                          periods_forecast,
                                          out_of_sample,xreg_out)
    new_table <- cbind(new_table,table)
  }
  get_table_min(new_table)
  return(new_table)
}

#generates model, then used autocompare
auto_forecast <- function(data,xreg_in,periods_forecast,
                           out_of_sample,xreg_out){
  ar0 <- c(1,0,1)
  new_table <- get_forecast_evaluation_exo(data,xreg_in,ar0,periods_forecast,
                              out_of_sample,xreg_out)
  i <- -1
  p <- 2
  q <- 2
  d <- 3#put back to 0
  while(i <= 5){
    new_p <- p + i
    j <- -1
    while(j <= 5){
      new_q <- q + j
      new_list_parameters <- c(new_p,d,new_q)
      table <- get_forecast_evaluation_exo(data,xreg_in,new_list_parameters,periods_forecast,
                                  out_of_sample,xreg_out)
      new_table <- cbind(new_table,table)
      j <- j + 1
    }
    i <- i + 1
  }
  
  get_table_min(new_table)
  return(new_table)
}


####### Tests
#unless stated other wise, these functions replace implementations in "point forecast evaluation_3.1.1.r"

####    Mincer-Zarnowitz (MZ) regression
# Check whether forecasts are unbiased and (strongly) efficient
# If we do not reject the null hypothesis, we have point forecasts that are unbiased
# The closer p_value is to 1, the more likely there isn't a significant difference between the value of alpha and zero as well as beta and 1
# It also tells us if the forecast is (strongly) efficient
# yhat should be pseudo out of sample forecast
# mz_test_simple only test for 1-step-ahead
mz_test_simple <- function(y,yhat){
  model<-lm(y~yhat)
  result <- linearHypothesis(model,c("(Intercept)=0", "yhat=1"))$`Pr(>F)`[2]
  return(result)
}

mz_test <- function(pasta_dataframe,pasta_forecast_matrix,length_out_sample,length_in_sample,horizon){
  forecast_length <- length_out_sample - horizon
  #forecast_length is due to remove first values when making higher h forecasts
  starting_forecast <- length_in_sample + 1
  #starting_forecast is due to forcasting the first value after length_in_sample
  ending_forecast <- length_out_sample + length_in_sample - horizon
  #ending_forecast is due to to the fact that we start after length_in_sample, but we remove the first values when making higher h forecasts
  modelMZ<-c()
  vector_result_MZ<-c()
  for (h in 1:horizon) {
    y <- pasta_dataframe[(starting_forecast+h):(ending_forecast+h)]
    yhat <- pasta_forecast_matrix[1:forecast_length,h]
    modelMZ[[h]]<- lm(y ~ yhat)
    
    vector_result_MZ[h]<-linearHypothesis(modelMZ[[h]], c("(Intercept)=0", "yhat=1"))$`Pr(>F)`[2]
  }
  return(vector_result_MZ)
}

#### Tau Test
## Direct way of testing forecast unbiasedness
## If tau=0, then unbiased
## The closer t_value is to 0, the more likely there isn't a significant difference between the value of tau and zero

# This function assumes the errors haven't been computed, and computes them
#It also assumes you are computing them for one horizon
# Might be a good idea to compute them beforehand for reporting
tau_test_simple <- function(y,yhat){
    
  error <- y - yhat
  model<-lm(error ~ 1)
  result <-coeftest(model, vcov = vcovHAC(model))
  
  tau_result<-matrix(ncol = 3)
  colnames(tau_result)<-c("coefficient","t-value", "p-value")
  tau_result[,1]<-result[1]#coefficient
  tau_result[,2]<-result[3]#t_vaue
  tau_result[,3]<-result[4]#p_value
  
  return(tau_result)
}

# Matrix of forecast errors
# CG should be changed to Neadly blabla
matrix_errors <- function(y,model_param,horizon){
  
  for1 <- function(x,h,parameters){
    forecast(Arima(x, order=c(parameters),optim.method = "Nelder-Mead"), h=h)}
  
  error_matrix <- tsCV(y,for1,h=horizon,parameters=model_param)
  
  return(error_matrix)
}
# Matrix of forecast errors, but model needs to be speficied by hand!!!!
# CG should be changed to Neadly blabla
matrix_errors_exo <- function(y,horizon,xreg){
  
  for1 <- function(x,h,parameters,xreg,newxreg){
    forecast(Arima(x, order=c(1,0,1),xreg=xreg), h=h,xreg=newxreg)}
  
  error_matrix <- tsCV(y,for1,h=horizon,xreg=xreg)
  
  return(error_matrix)
}

# This function assumes the forecast errors have been computed and are in matrix form
tau_test <- function(error_matrix,horizon){
  tau_result<-matrix(nrow = horizon,ncol = 3)
  colnames(tau_result)<-c("coefficient","t-value", "p-value")
  
  for(h in 1:horizon){
    model<-lm(error_matrix[,h] ~ 1)
    result <-coeftest(model, vcov = vcovHAC(model))
    
    tau_result[h,1]<-result[1]#coefficient
    tau_result[h,2]<-result[3]#t_vaue
    tau_result[h,3]<-result[4]#p_value
  }
  
  return(tau_result)
}

####### Spaghetti
### Reason behind the names: Spaghetti is a specific type of pasta
### All the functions with "spaghetti" only work on a specific code, when I tried
### to apply them to another scripts, even those provided by the courses, the
### functions did not work. 
### I then re wrote them in order to make them work everywhere. 
### So the name pasta because it is more general
pasta_dataframe<-function(Spaghettis_dataframe_initial,
                          length_out_sample, length_in_sample,horizon,model){
  
  len_initial_df <- length(Spaghettis_dataframe_initial[,2])
  rows_final <- len_initial_df+horizon #the data + the forecasts
  cols_final <- length_out_sample+2 #date, data and each step in forcast
  new_pasta_df <- data.frame(matrix(NA, nrow = rows_final, ncol = cols_final))
  new_pasta_df[1:len_initial_df,1] <- Spaghettis_dataframe_initial[,1]
  new_pasta_df[1:len_initial_df,2] <- Spaghettis_dataframe_initial[,2]
  for (c in 1:length_out_sample) {
    print(c)
    libor.fit <- Arima(Spaghettis_dataframe_initial[1:(length_in_sample + c), 2], order = c(model))
    new_pasta_df[(length_in_sample+1+c):(length_in_sample+c+horizon),2+c] <- forecast(libor.fit, h = horizon)$mean
    
    
  }
  new_pasta_df
}

pasta_dataframe_exo<-function(Spaghettis_dataframe_initial,
                              length_out_sample, length_in_sample,horizon,model,xreg){
  len_initial_df <- length(Spaghettis_dataframe_initial[,2])
  rows_final <- len_initial_df+horizon #the data + the forecasts
  cols_final <- length_out_sample+2 #date, data and each step in forcast
  new_pasta_df <- data.frame(matrix(NA, nrow = rows_final, ncol = cols_final))
  new_pasta_df[1:len_initial_df,1] <- Spaghettis_dataframe_initial[,1]
  new_pasta_df[1:len_initial_df,2] <- Spaghettis_dataframe_initial[,2]
  for (c in 1:length_out_sample) {
    #print(c)
    libor.fit <- Arima(Spaghettis_dataframe_initial[1:(length_in_sample + c), 2], order = c(model),xreg=xreg[1:(length_in_sample + c)])
    new_pasta_df[(length_in_sample+1+c):(length_in_sample+c+horizon),2+c] <- forecast(libor.fit, h = horizon,xreg=xreg[(length_in_sample+1+c):(length_in_sample+c+horizon)])$mean
    
    
  }
  new_pasta_df
}

pasta_plot<- function(Spaghetti_dataframe,
                      length_out_sample,  length_in_sample, horizon){
  
  #Spaghetti_dataframe= dataframe of the format of the above function Spaghetti_dataframe()
  
  plot(Spaghetti_dataframe[1:(length_in_sample+length_out_sample+ horizon+1), 2], type="l", xaxt="n", xlab="",ylab="" )
  for (i in 1:length_out_sample){
    
    lines(Spaghetti_dataframe[1:(length_in_sample+length_out_sample+ horizon +1),(2+i)], type = "l", col="lightblue")
    
  }
  
}

pasta_forecast_matrix<- function( spaghetti_dataframe, length_out_sample,
                                  length_in_sample, horizon){
  #*better practice to create "matrix_empty" inside function and output it.
  #*principle of functional programming, functions should not modify existing objects
  #matrix_empty= matrix with nrow= out sample and ncol = horizon
  #spaghetti dataframe= dataframe of format after function spaghetti_dataframe()
  matrix_empty <- matrix(nrow = length_out_sample, ncol = horizon)
  for (c in 1:length_out_sample){
    
    matrix_empty[c,1:horizon]<- spaghetti_dataframe[(length_in_sample+1+c):(length_in_sample+c+horizon),2+c]
  }
  matrix_empty
}

#### IC functions
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
  
  print("MIN BIC")
  print(min(matrix.BIC))
  print(matrix.BIC)
  
  print("MIN HQ")
  print(min(matrix.HQ))
  print(matrix.HQ)
  
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

