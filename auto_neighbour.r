library(forecast)

#compute model
get_ic <- function(data,p,d,q){
  model = Arima(data, order = c(p,d,q))		
  model_aic <- model$aicc
  model_bic <- model$bic
  names_ic <- c("aic","bic")
  list_ic <- c(model_aic,model_bic)
  name_col <- paste(c(p,d,q), collapse = "_")
  final_list <- as.data.frame(list_ic,row.names=names_ic)
  colnames(final_list) =name_col
  return(final_list)
}
get_neighbours <- function(data_table){
  print(data_table)
  print("Min AIC:")
  print(min(data_table[1,]))
  print("Min BIC:")
  print(min(data_table[2,]))
}

auto_neighbour <- function(data, p,d,q){
  if(p == 0){
    i <- 0
  }
  else{
    i <- -1
  }
  table <- get_ic(data,p,d,q)
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
      new_table <- get_ic(data,new_p,d,new_q)
      table <- cbind(table,new_table)
      j <- j + 1
    }
    i <- i + 1
  }

  get_neighbours(table)
  return(table)
}


auto_compare <- function(data,list_parameters){
  table <- get_ic(data)
}

library(MLmetrics)
get_forecast_evaluation <- function (data,parameter,periods,out_sample) {
  fit_basic1 <- Arima(data, order=c(parameter))
  forecast_1<-forecast(fit_basic1, h=periods)
  #lot(forecast_1)
  #plot(forecast_1$residuals)
  #acf(forecast_1$residuals)
  #pacf(forecast_1$residuals)
  summary(fit_basic1)
  #qqnorm(forecast_1$residuals)
  print("Model:")
  print(parameter)
  print("MAPE:")
  print(MAPE(forecast_1$mean, out_sample) * 100)
}


