library(fitdistrplus) # fitdist
library(gamlss.dist) # GG
library(zoo) # rollapply
library(stringr) # str_detect
library(magrittr) # %>%
library(TSstudio) # ts_split
library(forecast) # auto.arima, fourier, forecast
library(corrplot) # corrplot
library(vars) # VAR

load('traffic.RData')

TRAIN_RATIO <- 0.7
types <- c('UPLOAD', 'DOWNLOAD')
timestamps <- sort(unique(traffic_data$START_HOUR))
app_class <- as.vector(unique(traffic_data$APP_CLASS))
app_color <- c('chartreuse', 'coral2', 'cyan2', 'darkgoldenrod1', 'blue', 'forestgreen', 'darkviolet',
               'deeppink', 'mediumspringgreen', 'olivedrab2', 'red3', 'orangered', 'slateblue3', 'magenta',
               'brown', 'dimgrey')

'%*%' <- function(s1, s2)
{
  res <- c()
  for(j in s2)
    for(i in s1)
      res <- c(res, paste0(i,j))
  return(res)
}

### GG parameters and anomalies ###

reparamGG <- function(param)
{
  mu0 <- param[1]
  sigma <- param[2]
  nu <- param[3]
  r <- 1/(sigma*nu)^2
  return( list(r = as.numeric(r), gamma = as.numeric(nu), mu = as.numeric(r/mu0^nu)) )
}

detect_all_anomalies <- function(app, time, type, alpha=0.05)
{
  par(cex.main = 3, cex.lab = 2, cex.axis = 2)
  if(app == 'Full Traffic'){
    data <- traffic_data[traffic_data$START_HOUR %in% time, type] 
  }else{
    data <- traffic_data[traffic_data$START_HOUR %in% time & traffic_data$APP_CLASS == app, type] 
  }
  data <- data[data > 0]
  m <- length(data)
  date <- strsplit(time[1], ' ')[[1]][1]
  start_time <- (strsplit(time[1], ' ')[[1]][2] %>% strsplit(., ':'))[[1]][1]
  if(length(time) == 1){
    name <- paste(app, type, time) 
  }else{
    name <- paste(app, type, date, paste0(start_time, '-', sprintf('%02i', as.numeric(start_time)+4)))
  }
  fd <- fitdist(data, 'GG', method = 'mle', start = list(mu=1, sigma=2, nu=1), lower = c(0,0,-Inf))
  param <- reparamGG(fd$estimate)
  print(param$gamma)
  plot(data, ylab = ' ', xlab = 'Index', type = 'h', main = name)
  mtext('Value', side = 2, line = 2.85, cex = 2)
  for(i in 1:m){
    SR0 <- (m-1)*data[i]^param$gamma/sum(data[-i]^param$gamma)
    if(param$gamma > 0){
      pval <- 1 - pf(SR0, df1 = param$r, df2 = (m-1)*param$r)
    }else{
      pval <- 1 - pf(1/SR0, df1 = (m-1)*param$r, df2 = param$r)
    }
    if(pval < alpha){
      points(i, data[i], col = 'red', pch = 18, cex = 2.5)
    }
  }
  text(x = m, y = 10/11*max(data), labels = bquote(gamma == .(round(param$gamma,3))), pos = 2, cex = 4)
  par(cex.main = 1, cex.lab = 1, cex.axis = 1)
}

plot_ggamma <- function(app, time, type)
{
  par(cex.main = 3, cex.lab = 1.5, cex.axis = 2)
  if(any(str_detect(time, ' '))){
    time_cond <- traffic_data$START_HOUR %in% time
  }else{
    time_cond <- str_detect(traffic_data$START_HOUR, time)
  }
  if(app == 'Full Traffic'){
    data <- traffic_data[time_cond, type] 
  }else{
    data <- traffic_data[time_cond & traffic_data$APP_CLASS == app, type]  
  }
  data <- data[data > 0]
  nobs <- length(data)
  m <- max(data)
  q <- quantile(data, 0.75) + 1.75*IQR(data)
  fd <- fitdist(data, 'GG', method = 'mle', start = list(mu=1, sigma=2, nu=1), lower = c(0,0,-Inf))
  x <- seq(1e-32, 1.1*q, by = q/10000)
  y <- sapply(x, function(z) dGG(z, mu = fd$estimate[1], sigma = fd$estimate[2], nu = fd$estimate[3]))
  date <- strsplit(time[1], ' ')[[1]][1]
  start_time <- (strsplit(time[1], ' ')[[1]][2] %>% strsplit(., ':'))[[1]][1]
  if(length(time) == 1){
    name <- paste(app, type, time) 
  }else{
    name <- paste(app, type, date, paste0(start_time, '-', sprintf('%02i', as.numeric(start_time)+4)))
  }
  br <- 'FD'
  if(app %in% c('VoIP', 'Games', 'Others')){
    br <- seq(0, m + q/16, by = q/16)
  }
  if(app == 'Full Traffic'){
    #br <- seq(from = 0, to = m, length.out = 500)
    #br <- 'Scott'
    #br <- 'FD'
    #q_arg <- seq(0,1,length.out=len)
    #br <- qGG(q_arg, mu = fd$estimate[1], sigma = fd$estimate[2], nu = fd$estimate[3])
    #br <- quantile(data, q_arg)
    #br <- ceiling(nobs^(1/3)*m/(3.5*sd(data)))
    #br <- (5*log(nobs,10) - 5) %>% round()
    br <- sqrt(nobs) %>% round()
  }
  if(app == 'Streaming Applications'){
    br <- seq(0, m + q/64, by = q/64)
  }
  maxdens <- max(hist(data, breaks = br, plot = FALSE)$density)
  if(is.nan(y[1]) | y[1] < 0.5){
    ymax <- max(c(y[!is.nan(y)], maxdens))
  }else{
    ymax <- maxdens
  }
  hm <- hist(data, freq = FALSE, breaks = br, xlim = c(0,q/2),  ylim = c(0,ymax), xlab = 'Value', main = name)
  lines(x, y, col = 'blue', lwd = 2)
  br_cdf <- do.call(paste0('p', fd$distname), c(list(hm$breaks), lapply(fd$estimate, function(x) x)))
  probs <- rollapply(br_cdf, 2, function(x) x[2]-x[1])
  chi <- chisq.test(hm$counts, p = probs, rescale.p = TRUE)#, simulate.p.value = TRUE, B=1000) # Monte Carlo simulations
  text(x = q/6, y = 4*ymax/11, labels = paste0('Mu = ', round(fd$estimate[1], 4), ', Sigma = ', round(fd$estimate[2], 4), ', Nu = ', round(fd$estimate[3], 4)), pos = 4, cex = 2)
  text(x = q/6, y = 4.5*ymax/11, labels = paste(nobs, 'Observations'), pos = 4, cex = 2)
  text(x = q/6, y = 5*ymax/11, labels = paste('Chi-squared test p-value ', round(chi$p.value,4)), pos = 4, cex = 2)
  par(cex.main = 1, cex.lab = 1, cex.axis = 1)
  return(c(fd$estimate, chi$p.value))
}

### Forecasting ###

## 1 hour ##
l <- 24 # frequency
s <- c(10,2) # start time

## 4 hour ##
l <- 6 # frequency
s <- c(10,1) # start time

cor.mtest <- function(mat, ...) 
{
  mat <- as.matrix(mat)
  n <- ncol(mat)
  p.mat<- matrix(NA, n, n)
  diag(p.mat) <- 0
  for (i in 1:(n-1)) {
    for (j in (i+1):n) {
      tmp <- cor.test(mat[ ,i], mat[ ,j], ...)
      p.mat[i,j] <- p.mat[j,i] <- tmp$p.value
    }
  }
  colnames(p.mat) <- rownames(p.mat) <- colnames(mat)
  return(p.mat)
}

p.mat <- cor.mtest(mat) # mat - time series (columns)  
col_panel <- colorRampPalette(c("#BB4444", "#EE9988", "#FFFFFF", "#77AADD", "#4477AA"))
mat %>% cor() %>% corrplot(., method = 'color', col=col_panel(200),  
                             type='upper', 
                             addCoef.col = 'black', 
                             tl.col = 'black', tl.srt=45, 
                             p.mat = p.mat, sig.level = 0.01, insig = 'pch', 
                             pch.col = 'red', diag=FALSE,
                             number.cex = 1.9, tl.cex = 2, cl.cex = 2, pch.cex = 7)

rmse <- function(x, y)
{
  if(length(x) != length(y)){
    print("Vectors are not the same length!")
    return (NA)
  }
  return ( sqrt( mean( (x - y)^2 ) ) )
}

calc_rmse_train <- function(ts, kfourier)
{
  ntest <- round(length(ts)*(1 - TRAIN_RATIO))
  split <- ts_split(ts, sample.out = ntest)
  train <- split$train
  m <- auto.arima(train, trace = TRUE, approximation = TRUE, stepwise = TRUE, xreg = fourier(train, K = kfourier))
  return(rmse(train, m$fitted))
}

calc_rmse_var <- function(ts, p)
{
  l <- frequency(ts)
  ntest <- round(dim(ts)[1]*(1-TRAIN_RATIO))
  split <- ts_split(ts, sample.out = ntest)
  train <- split$train
  test <- split$test
  trmin <- sapply(1:dim(ts)[2], function(i) min(train[ ,i]))
  trmax <- sapply(1:dim(ts)[2], function(i) max(train[ ,i]))
  for(i in 1:dim(ts)[2]){
    train[ ,i] <- (train[ ,i] - trmin[i])/(trmax[i] - trmin[i])
    test[ ,i] <- (test[ ,i] - trmin[i])/(trmax[i] - trmin[i])
  }
  m <- VAR(train, p = p, season = l, exogen = NULL, type = 'const')
  f <- forecast(m, h = ntest, level = 95)
  rmse_vec <- rep(NA, dim(ts)[2])
  for(i in 1:dim(ts)[2]){
    rmse_vec[i] <- rmse(test[ ,i], f$forecast[[i]]$mean)
  }
  return(rmse_vec)
}

plot_var <- function(ts, p)
{
  l <- frequency(ts)
  ntest <- round(dim(ts)[1]*(1-TRAIN_RATIO))
  split <- ts_split(ts, sample.out = ntest)
  train <- split$train
  test <- split$test
  start_time <- Sys.time()
  m <- VAR(train, p = p, season = l, exogen = NULL, type = 'const')
  f <- forecast(m, h = ntest, level = 95)
  dur <- as.numeric(difftime(Sys.time(), start_time, units = 'secs'))
  rmse_vec <- rep(NA, dim(ts)[2])
  par(cex.main = 3, cex.lab = 1.5, cex.axis = 1.5)
  for(i in 1:dim(ts)[2]){
    plot(f$forecast[[i]], lwd = 2, fcol = 'burlywood1', shadecols = 'burlywood1', type = 'b', ylim = c(min(c(ts[ ,i], f$forecast[[i]]$lower)),max(c(ts[ ,i], f$forecast[[i]]$upper))), xlab = 'Time')
    lines(test[ ,i], col = 'black', lwd = 2, type = 'b')
    lines(f$forecast[[i]]$mean, col = 'orangered', lwd = 2, type = 'b')
    abline(v = end(train)[1] + (end(train)[2]-1)/frequency(train), lty = 2, lwd = 2, col = 'orangered')
    text(x = 11, y = min(ts[ ,i]), labels = paste0('RMSE=', round(rmse(test[ ,i], f$forecast[[i]]$mean), 5)), pos = 4, cex = 2)
    abline(v = end(train)[1] + (end(train)[2]-1)/frequency(train), lty = 2, lwd = 2, col = 'orangered')
    rmse_vec[i] <- rmse(test[ ,i], f$forecast[[i]]$mean)
  }
  par(cex.main = 1, cex.lab = 1, cex.axis = 1)
  return(dur)
}

calc_rmse_auto_arima <- function(ts, kfourier=0){
  ntest <- round(length(ts)*(1-TRAIN_RATIO))
  split <- ts_split(ts, sample.out = ntest)
  train <- split$train
  test <- split$test
  trmin <- min(train)
  trspread <- max(train) - min(train)
  norm_train <- (train - trmin)/trspread
  norm_test <- (test - trmin)/trspread
  if(kfourier == 0){
    m <- auto.arima(norm_train, trace = TRUE, approximation = FALSE, stepwise = TRUE)
    f <- forecast(m, h = ntest, level = 95)
  }else{
    m <- auto.arima(norm_train, trace = TRUE, approximation = FALSE, stepwise = TRUE, xreg = fourier(norm_train, K=kfourier))
    f <- forecast(m, xreg = fourier(norm_train, K=kfourier, h=ntest), level = 95)
  }
  return(rmse(norm_test, f$mean))
}

plot_auto_arima <- function(ts, title, kfourier=0)
{
  par(cex.main = 3, cex.lab = 1.5, cex.axis = 1.5)
  l <- frequency(ts)
  ntest <- round(length(ts)*(1-TRAIN_RATIO))
  plot(ts, type = 'b', lwd = 2, ylab = '', main = title)
  abline(reg = lm(ts ~ time(ts)), col = 'darkblue', lw = 2)
  split <- ts_split(ts, sample.out = ntest)
  train <- split$train
  test <- split$test
  start_time <- Sys.time()
  if(kfourier == 0){
    m <- auto.arima(train, trace = TRUE, approximation = FALSE, stepwise = TRUE)
    f <- forecast(m, h = ntest, level = 95)
    dur <- difftime(Sys.time(), start_time, units = 'secs')
  }else{
    m <- auto.arima(train, trace = TRUE, approximation = FALSE, stepwise = TRUE, xreg = fourier(train, K=kfourier))
    f <- forecast(m, xreg = fourier(train, K=kfourier, h=ntest), level = 95)
    dur <- difftime(Sys.time(), start_time, units = 'secs')
  }
  plot(f, lwd = 2, fcol = 'burlywood1', shadecols = 'burlywood1', type = 'b', ylim = c(min(c(ts, f$lower)), max(c(ts, f$upper))), xlab = 'Time')
  lines(test, col = 'black', lwd = 2, type = 'b')
  lines(f$mean, col = 'orangered', lwd = 2, type = 'b')
  abline(v = end(train)[1] + (end(train)[2]-1)/frequency(train), lty = 2, lwd = 2, col = 'orangered')
  text(x = 11, y = min(ts), labels = paste0('RMSE=', round(rmse(test, f$mean), 5)), pos = 4, cex = 2)
  par(cex.main = 1, cex.lab = 1, cex.axis = 1)
  return(as.numeric(dur))
}

### Users Sum Avg ###

users_u <- sapply(timestamps, function(t) length(unique(traffic_data[traffic_data$START_HOUR == t & traffic_data$UPLOAD > 0 ,'MASKED_MSISDN'])))
users_d <- sapply(timestamps, function(t) length(unique(traffic_data[traffic_data$START_HOUR == t & traffic_data$DOWNLOAD > 0 ,'MASKED_MSISDN']))) 
users <- cbind(users_u, users_d)

plot_users_hist <- function(users, type)
{
  par(mfrow = c(2,1), cex.main = 3, cex.lab = 2, cex.axis = 2)
  traffic_types <- c('UPLOAD', 'DOWNLOAD')
  tp <- match(type, traffic_types)
  nobs <- dim(users)[1]
  n <- which(row.names(users) == '2018-02-22 13:00:00')
  ints <- list(1:n, (n+1):nobs)
  for(i in ints){
    h <- hist(users[i,tp], freq = FALSE, xlab = 'Number of Unique Users', ylab = '',
              main = paste(row.names(users)[i][1], '-', row.names(users)[i][length(i)], type))
  }
  ks_pval <- ks.test(users[ints[[1]],tp], users[ints[[2]],tp])$p.value
  par(mfrow = c(1,1))
  grid <- seq(0.8*min(users[ ,tp]), 1.1*max(users[ ,tp]), by = (max(users[ ,tp])-min(users[ ,tp]))/10000)
  plot(ecdf(users[ints[[1]],tp]), col = 'deepskyblue', lwd = 2, cex = 1.3, 
       xlab = 'Number of Unique Users', ylab = '', main = paste('CDF', type))
  lines(ecdf(users[ints[[2]],tp]), col = 'orange', lwd = 2, cex = 1.3)
  text(x = min(users[ints[[1]],tp]), y = 0.6, labels = paste('K-S test p-value', round(ks_pval,4)), pos = 4, cex = 2)
  legend('bottomright', c('Train ecdf', 'Test ecdf'), col = c('deepskyblue', 'orange'), lwd = 3, cex = 1.5, text.width = 2000)
  par(cex.main = 1, cex.lab = 1, cex.axis = 1)
  return(ks_pval)
}

calc_hourly_users <- function(data, type) 
{ sapply(timestamps, function(t) length(unique(data[data$START_HOUR == t & data[ ,type] > 0, 'MASKED_MSISDN']))) }

calc_hourly_users_all <- function(data) 
{ sapply(timestamps, function(t) length(unique(data[data$START_HOUR == t, 'MASKED_MSISDN']))) }

calc_hourly_sum <- function(data, type)
{ sapply(timestamps, function(t) sum(data[data$START_HOUR == t, type])) }

collect_hourly_info <- function(app)
{
  cur_matr <- matrix(NA, nrow = length(timestamps), ncol = 9)
  row.names(cur_matr) <- timestamps
  colnames(cur_matr) <- c('u_users', 'u_sum', 'u_avg', 'd_users', 'd_sum', 'd_avg', 'all_users', 'all_sum', 'all_avg')
  data <- traffic_data[traffic_data$APP_CLASS == app, c('START_HOUR', 'MASKED_MSISDN', 'UPLOAD', 'DOWNLOAD')]
  cur_matr[ ,'u_users'] <- calc_hourly_users(data, 'UPLOAD')
  cur_matr[ ,'d_users'] <- calc_hourly_users(data, 'DOWNLOAD')
  cur_matr[ ,'all_users'] <- calc_hourly_users_all(data)
  cur_matr[ ,'u_sum'] <- calc_hourly_sum(data, 'UPLOAD')
  cur_matr[ ,'d_sum'] <- calc_hourly_sum(data, 'DOWNLOAD')
  cur_matr[ ,'all_sum'] <- cur_matr[ ,'u_sum'] + cur_matr[ ,'d_sum']
  cur_matr[ ,'u_avg'] <- cur_matr[ ,'u_sum']/cur_matr[ ,'u_users']
  cur_matr[ ,'d_avg'] <- cur_matr[ ,'d_sum']/cur_matr[ ,'d_users']
  cur_matr[ ,'all_avg'] <- cur_matr[ ,'all_sum']/cur_matr[ ,'all_users']
  return(cur_matr)
}

app_hour_param <- list()
for(app in app_class){
  app_hour_param[[app]] <- collect_hourly_info(app)
}

col_alias <- c('Number of Unique Users ', 'Total ', 'Average per User ') %*% c('UPLOAD Traffic', 'DOWNLOAD Traffic', 'ALL Traffic')
names(col_alias) <- colnames(app_hour_param[[1]])

plot_hourly_info <- function(colname)
{
  par(cex.main = 3, cex.lab = 2, cex.axis = 1.5)
  nobs <- dim(app_hour_param[[1]])[1]
  n <- which(rownames(app_hour_param[[1]]) == '2018-02-22 13:00:00')
  ints <- list(1:n, (n+1):nobs)
  matr <- matrix(NA, nrow = dim(app_hour_param[[1]])[1], ncol = length(app_class))
  rownames(matr) <- rownames(app_hour_param[[1]])
  colnames(matr) <- app_class
  for(app in app_class)
    matr[ ,app] <- app_hour_param[[app]][ ,colname]
  matr[is.nan(matr)] <- 0
  plot(1:nobs, rep(NA, nobs), xaxt='n', xlab = '', ylim = c(min(matr), max(matr)), ylab = '', main = col_alias[colname])
  axis(1, at=1:nobs, labels = rownames(app_hour_param[[1]]))
  abline(v=n, col = 'darkgray', lty = 2, lwd = 2)
  abline(v=n+1, col = 'darkgray', lty = 2, lwd = 2)
  for(j in 1:length(app_class)){
    for(i in ints){
      lines(i, matr[i,app_class[j]], col = app_color[j], lwd = 2, type = 'l')
    }
  }
  par(cex.main = 1, cex.lab = 1, cex.axis = 1)
}

app_params <- colnames(app_hour_param[[1]])
app_hour_param_by_params <- list()
for(item in app_params){
  cur_matr <- matrix(NA, nrow = length(timestamps), ncol = length(app_class))
  rownames(cur_matr) <- timestamps
  colnames(cur_matr) <- app_class
  for(app in app_class){
    cur_matr[ ,app] <- app_hour_param[[app]][ ,item]
  }
  app_hour_param_by_params[[item]] <- cur_matr
}

mark <- 301
test_length <- length(timestamps) - mark
nahead <- 10 + 6*24 + 4*24

calc_rmse_var_after_break <- function(ts, p, nahead, ntest)
{
  l <- frequency(ts)
  split <- ts_split(ts, sample.out = ntest)
  train <- split$train
  test <- split$test
  trmin <- sapply(1:dim(ts)[2], function(i) min(train[ ,i]))
  trmax <- sapply(1:dim(ts)[2], function(i) max(train[ ,i]))
  for(i in 1:dim(ts)[2]){
    train[ ,i] <- (train[ ,i] - trmin[i])/(trmax[i] - trmin[i])
    test[ ,i] <- (test[ ,i] - trmin[i])/(trmax[i] - trmin[i])
  }
  m <- VAR(train, p = p, season = l, exogen = NULL, type = 'const')
  f <- forecast(m, h = nahead, level = 95)
  rmse_vec <- rep(NA, dim(ts)[2])
  len <- length(f$forecast[[1]]$mean)
  for(i in 1:dim(ts)[2]){
    rmse_vec[i] <- rmse(test[ ,i], f$forecast[[i]]$mean[(len-ntest+1):len])
  }
  return(rmse_vec)
} 
