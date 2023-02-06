library(quantmod)
library(xts)
library(zoo,warn.conflicts=FALSE)
library(lubridate,warn.conflicts=FALSE)
library(mgcv,warn.conflicts=FALSE)
library(rugarch,warn.conflicts=FALSE)
## load the data
getSymbols("^DJI",src="yahoo")
mydates = row.names(as.data.frame(DJI))
mydates=as.Date(mydates,"%Y-%m-%d")
candleChart(DJI,multi.col=TRUE,theme="white")
dj <- Cl(DJI)
dj=xts(dj,mydates)

## plots
plot(dj, main='DJI Closed Price')
plot(diff(dj), main='DJI 1 lag Differenced Price')

## non-parametric trend using splines regression
time.pts = c(1:length(dj))
time.pts = c(time.pts - min(time.pts))/max(time.pts)
dj.gam = gam(dj~s(time.pts))
dj.fit = xts(fitted(dj.gam),mydates)
plot(dj, ylab='Price', main='Dow Jones Industrial Average')
lines(dj.fit, lwd=3, col='orange')

## alternative method, using 'ts'
dj <- ts(DJI$DJI.Close,start = c(2007,1,1),freq=252)
dj.gam = gam(dj~s(time.pts))
dj.fit = ts(fitted(dj.gam),start = c(2007,1,1),freq=252)
ts.plot(dj,ylab='Price', main='Dow Jones Industrial Average')
lines(dj.fit, lwd=2, col='red')
grid()


## return closed price (following analysis will use the return price)
dj.return <- diff(dj) / dj[-length(dj)]
plot(dj.return, main='Daily Return Price')
grid()
acf(dj.return, main='ACF of Daily Return Price')
pacf(dj.return, main='PACF of Daily Return Price')

## ARIMA
  ## Split data into training and testing sets
  n <- length(dj.return)
  nfit <- n-8
  dj.train <- dj.return[1:nfit]
  dj.test <- dj.return[(nfit+1):n]

  ## Fit ARIMA = (p,1,q)
  norder <- 8
  p <- c(1:norder)-1
  q <- c(1:norder)-1
  dj.aic1 <- matrix(0,norder,norder)
  for(i in 1:norder){
    for(j in 1:norder){
      modij = arima(dj.train,order = c(p[i],1,q[j]), seasonal = list(order = c(0,0,0)),
                    method='ML',optim.control = list(maxit = 1000))
      dj.aic1[i,j] = modij$aic-2*(p[i]+q[j]+1)+2*(p[i]+q[j]+1)*nfit/(nfit-p[i]-q[j]-2)
    }
  }

  ## Fit ARIMA = (p,2,q)
  norder <- 8
  p <- c(1:norder)-1
  q <- c(1:norder)-1
  dj.aic2 <- matrix(0,norder,norder)
  for(i in 1:norder){
    for(j in 1:norder){
      modij = arima(dj.train,order = c(p[i],2,q[j]), method='ML')
      dj.aic2[i,j] = modij$aic-2*(p[i]+q[j]+1)+2*(p[i]+q[j]+1)*nfit/(nfit-p[i]-q[j]-2)
    }
  }
  
  ## Final ARIMA model
  if (min(dj.aic1)==min(dj.aic1,dj.aic2)){ 
    dj.aic <- dj.aic1
    dorder <- 1
  } else{
    dj.aic <- dj.aic2
    dorder <- 2
  }
  
  porder <- as.numeric(which(dj.aic==min(dj.aic),arr.ind = TRUE)[1,1]-1)
  qorder <- as.numeric(which(dj.aic==min(dj.aic),arr.ind = TRUE)[1,2]-1)
  
  ## Final Model Fit
  final.dj <- arima(dj.train, order = c(porder,dorder,qorder), method='ML')
  summary(final.dj)
  final.dj$coef
  final.dj$var.coef
  final.dj$aic
  cbind(porder, dorder, qorder)
  
  # Residuals
  acf(final.dj$residuals)
  pacf(final.dj$residuals)
  Box.test(final.dj$residuals, lag = (porder+qorder), 
           type = 'Ljung', fitdf = (porder+qorder))
  Box.test((final.dj$residuals)^2, lag = (porder+qorder), 
           type = 'Ljung', fitdf = (porder+qorder))

## ARCH


## GARCH
  #ARMA order
  #Fit ARMA on return price
  final.aic=Inf
  final.order=c(0,0,0)
  for (p in 0:9) for (d in 0:0) for (q in 0:9) {
    current.aic=AIC(arima(dj.train,order=c(p, d, q), method='ML'))
    if(current.aic<final.aic) {
      final.aic=current.aic
      final.order=c(p,d,q)
      final.arima=arima(dj.train, order=final.order)
    }
  }
  ## Find GARCH Order given ARMA order identified before
  ## ugrach from rugarch libary
  library(rugarch)
  #Initial GARCH Order
  #ARIMA-GARCH: Select GARCH order
  test_modelAGG <- function(m,n){
    spec = ugarchspec(variance.model=list(garchOrder=c(m,n)),
                      mean.model=list(armaOrder=c(9,9),
                                      include.mean=T),
                      distribution.model="std")
    fit = ugarchfit(spec, dj.train, solver = 'hybrid')
    current.bic = infocriteria(fit)[2]
    df = data.frame(m,n,current.bic)
    names(df) <- c("m","n","BIC")
    print(paste(m,n,current.bic,sep=" "))
    return(df)
  }
  
  ordersAGG = data.frame(Inf,Inf,Inf)
  names(ordersAGG) <- c("m","n","BIC")
  
  for (m in 0:4){
    for (n in 0:4){
      possibleError <- tryCatch(
        ordersAGG<-rbind(ordersAGG,test_modelAGG(m,n)),
        error=function(e) e
      )
      if(inherits(possibleError, "error")) next
    }
  }
  ordersAGG <- ordersAGG[order(-ordersAGG$BIC),]
  tail(ordersAGG)
  # 2,1
  
  #ARMA update
  test_modelAGA <- function(p,q){
    spec = ugarchspec(variance.model=list(garchOrder=c(2,1)),
                      mean.model=list(armaOrder=c(p,q),
                                      include.mean=T),
                      distribution.model="std")
    fit = ugarchfit(spec, dj.train, solver = 'hybrid')
    current.bic = infocriteria(fit)[2]
    df = data.frame(p,q,current.bic)
    names(df) <- c("p","q","BIC")
    print(paste(p,q,current.bic,sep=" "))
    return(df)
  }
  
  ordersAGA = data.frame(Inf,Inf,Inf)
  names(ordersAGA) <- c("p","q","BIC")
  for (p in 0:9){
    for (q in 0:9){
      possibleError <- tryCatch(
        ordersAGA<-rbind(ordersAGA,test_modelAGA(p,q)),
        error=function(e) e
      )
      if(inherits(possibleError, "error")) next
    }
  }
  ordersAGA <- ordersAGA[order(-ordersAGA$BIC),]
  tail(ordersAGA)
  #0,1
  
  #GARCH update
  test_modelAGG <- function(m,n){
    spec = ugarchspec(variance.model=list(garchOrder=c(m,n)),
                      mean.model=list(armaOrder=c(1,0),
                                      include.mean=T), distribution.model="std")
    fit = ugarchfit(spec, dj.train, solver = 'hybrid')
    current.bic = infocriteria(fit)[2]
    df = data.frame(m,n,current.bic)
    names(df) <- c("m","n","BIC")
    print(paste(m,n,current.bic,sep=" "))
    return(df)
  }
  
  ordersAGG = data.frame(Inf,Inf,Inf)
  names(ordersAGG) <- c("m","n","BIC")
  
  for (m in 0:2){
    for (n in 0:2){
      possibleError <- tryCatch(
        ordersAGG<-rbind(ordersAGG,test_modelAGG(m,n)),
        error=function(e) e
      )
      if(inherits(possibleError, "error")) next
    }
  }
  ordersAGG <- ordersAGG[order(-ordersAGG$BIC),]
  tail(ordersAGG)
  
  spec.dj = ugarchspec(variance.model=list(garchOrder=c(2,1)),
                         mean.model=list(armaOrder=c(1, 0),
                                         include.mean=T), distribution.model="std")
  final.model = ugarchfit(spec.dj, dj.train, solver = 'hybrid')
  final.model
  