library(readxl)
library(forecast)
library(stats)
library(lmtest)

# Limpiar ambiente
rm(list=ls())

sifmi_sect <- read_excel("Analisis sifmi dic 2021.xlsx", sheet="VarxSector")
sifmi_sect <- sifmi_sect[c(1:23),]
sifmi_sect <- as.data.frame(sifmi_sect[c(1:23),c(2:length(sifmi_sect))])
sifmi_sect <- abs(sifmi_sect)
sifmi_sect <- t(sifmi_sect)
colnames(sifmi_sect)<- c("Total","11","21","22","23","31","32","33","43","46","48","49","51","52","53","54","55",
                         "56","61","62","71","72","81")
sifmi_sect <- ts(sifmi_sect,start=2016,frequency = 12)

# Revisar que las columnas tengan nombre

#### ARIMA
# Inicializar variables del loop
Selec_Arima <- data.frame()
Pronostico <- data.frame()
Pronostico_ets <- data.frame()

# Función que busca auto arima sin especificar un criterio
estimacion3<-function(sector){
  tryCatch({
    modelo2 <- auto.arima(log(sifmi_sect[,sector]), max.d=1, max.Q = 0, max.P = 0, max.D = 0, stepwise = FALSE, method="ML")
    #summary(modelo2)
    valores_iniciales <- coeftest(modelo2)
    valores <- valores_iniciales
    
    coef_sign <- as.integer(valores_iniciales[,4]<0.05)
    terminos_arma <- replace(coef_sign, coef_sign==1,NA)
    
    modelo1 <- arima(log(sifmi_sect[,sector]),
                     order = c(modelo2[["arma"]][1],modelo2[["arma"]][6],modelo2[["arma"]][2]), 
                     method="ML", transform.pars = FALSE, fixed = terminos_arma)
    valores <- coeftest(modelo1)
    
    # Estacionaridad mediante raices del polinomio AR
    raiz <- abs(polyroot(c(1,-modelo1$model$phi)))
    
    # No autocorrelacion del error mediante Ljung-Box (p.value > 0.05)
    LB <- checkresiduals(modelo1, plot= FALSE)
    
    if(nrow(valores)==1){
      valores <- data.frame(Terminos_ARMA = row.names(valores), t(valores[1:nrow(valores),1:ncol(valores)]), Metodo=3, 
                            LBest= LB$statistic, LBvalorp= LB$p.value)
    } else{
      valores <- data.frame(Terminos_ARMA = row.names(valores), valores[1:nrow(valores),1:ncol(valores)], Metodo = 3, 
                 LBest= LB$statistic, LBvalorp=LB$p.value)
    }
    
    #valores$Raices <- list(raiz)
    
    fc <- forecast(modelo1,h=1)
    
    return(list(Arima=valores,forecast=fc))
    
  },error=function(e){cat("ERROR en función 3:", colnames(sifmi_sect)[sector], "\n")
    valores=data.frame(Terminos_ARMA = row.names(valores_iniciales), valores_iniciales[1:nrow(valores_iniciales),1:ncol(valores_iniciales)], Metodo=4, 
               LBest= NA, LBvalorp= NA)
    return(list(Arima=valores, forecast=data.frame(Point.Forecast=NA,Lo.80=NA,Hi.80=NA,Lo.95=NA,Hi.95=NA)))
  })
}

# Función que remueve todos los terminos ar con probabilidad mayor a 0.05 en una sola acción
estimacion2<-function(sector){
  tryCatch({
    modelo2 <- auto.arima(log(sifmi_sect[,sector]), max.d=1, max.Q = 0, max.P = 0, max.D = 0, ic="bic", stepwise = FALSE, method="ML")
    #summary(modelo2)
    valores_iniciales <- lmtest::coeftest(modelo2)
    valores <- valores_iniciales
    
    coef_sign <- as.integer(valores_iniciales[,4]<0.05)
    terminos_arma <- replace(coef_sign, coef_sign==1,NA)
    
    modelo1 <- arima(log(sifmi_sect[,sector]),
                     order = c(length(modelo2$model$phi),modelo2$model$Delta,length(modelo2$model$theta)), 
                     method="ML", transform.pars = FALSE, fixed = terminos_arma)
    valores <- lmtest::coeftest(modelo1)
    #print(valores)
    
    # Estacionaridad mediante raices del polinomio AR
    raiz <- abs(polyroot(c(1,-modelo1$model$phi)))
    
    # No autocorrelacion del error mediante Ljung-Box (p.value > 0.05)
    LB <- forecast::checkresiduals(modelo1, plot= FALSE)
    
    if(nrow(valores)==1){
      valores <- data.frame(Terminos_ARMA = row.names(valores), t(valores[1:nrow(valores),1:ncol(valores)]), Metodo=2, 
                            LBest= LB$statistic, LBvalorp= LB$p.value)
    } else{
      valores <- data.frame(Terminos_ARMA = row.names(valores), valores[1:nrow(valores),1:ncol(valores)], Metodo = 2, 
                            LBest= LB$statistic, LBvalorp=LB$p.value)
    }
    #valores$Raices <- list(raiz)
    
    fc <- forecast(modelo1,h=1)
    
    return(list(Arima=valores,forecast=fc))
    
  },error=function(e){cat("ERROR en función 2:", colnames(sifmi_sect)[sector], "\n")
    return(estimacion3(sector))
  })
}

# Función que remueve terminos ar con probabilidades mayores a 0.05 de uno por uno
estimacion<-function(sector){
  tryCatch({
    modelo2 <- auto.arima(log(sifmi_sect[,sector]), max.d=1, max.P=0, max.Q = 0, max.D = 0, ic="bic", stepwise = FALSE, method="ML")
    #summary(modelo2)
    valores_iniciales <- coeftest(modelo2)
    valores <- valores_iniciales
    
    # Estacionaridad mediante raices del polinomio AR
    raiz <- abs(polyroot(c(1,-modelo2$model$phi)))
    
    # No autocorrelacion del error mediante Ljung-Box (p.value > 0.05)
    LB <- forecast::checkresiduals(modelo2, plot= FALSE)
    
    fc <- forecast(modelo2,h=1)
    
    x <- 0
    
    while (!max(valores[,4])<0.05) {
      coef_sign <- as.integer(valores_iniciales[,4]<sort(valores_iniciales[,4],partial=length(valores_iniciales[,4])-x)[length(valores_iniciales[,4])-x])
      terminos_arma <- replace(coef_sign, coef_sign==1,NA)
      
      modelo1 <- arima(log(sifmi_sect[,sector]), 
                       order = c(length(modelo2$model$phi),modelo2$model$Delta,length(modelo2$model$theta)), 
                       method="ML", transform.pars = FALSE, fixed = terminos_arma)
      valores <- coeftest(modelo1)
      #print(valores)
      
      # Realiza las mismas pruebas si entra en el loop
      raiz <- abs(polyroot(c(1,-modelo1$model$phi)))
      LB <- checkresiduals(modelo1, plot= FALSE)
      
      fc <- forecast(modelo1,h=1)
      
      x <- x+1
    }
    
    if(nrow(valores)==1){
      valores <- data.frame(Terminos_ARMA = row.names(valores), t(valores[1:nrow(valores),1:ncol(valores)]), Metodo=1, 
                            LBest= LB$statistic, LBvalorp= LB$p.value)
    } else{
      valores <- data.frame(Terminos_ARMA = row.names(valores), valores[1:nrow(valores),1:ncol(valores)], Metodo = 1, 
                            LBest= LB$statistic, LBvalorp=LB$p.value)
    }
    #valores$Raices <- list(raiz)
    
    return(list(Arima=valores,forecast=fc))
    
  },error=function(e){cat("ERROR en función 1:", colnames(sifmi_sect)[sector], "\n")
    return(estimacion2(sector))
  })
}

for (sector in c(1:23)) {
  resultado<-estimacion(sector)
  Selec_Arima <- rbind(Selec_Arima, data.frame(resultado$Arima,Sector=colnames(sifmi_sect)[sector]))
  Pronostico <- rbind(Pronostico, data.frame(resultado$forecast, Sector=colnames(sifmi_sect)[sector]))
}

Pronostico <- data.frame(exp(Pronostico[1:5]),row.names = colnames(sifmi_sect)[1:23])
sum(Pronostico$Point.Forecast[2:23],na.rm = TRUE)

#### ETS
for(sector in c(1:23)){
  fit_ets<- ets(sifmi_sect[,sector], model = "ZZZ",ic="bic", opt.crit="lik", allow.multiplicative.trend = TRUE)
  fcast_ets <- forecast(fit_ets, h=1)
  Pronostico_ets <- rbind(Pronostico_ets, data.frame(sector=colnames(sifmi_sect)[sector], fcast_ets))
}

write.csv(Pronostico,file="2111_pronost.csv", row.names=TRUE)
write.csv(Pronostico_ets,file="2111_pronostets.csv", row.names=FALSE)
