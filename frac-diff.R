setwd("/Users/jvsingh/work/fredmd")

#rm(list = ls())
library(xts)
library(reshape2)

df   <- read.csv(file = "2019-04.txt", header = FALSE, sep = ",", skip = 2, stringsAsFactors = FALSE)
temp <- read.csv(file = "2019-04.txt", header = TRUE, sep = ",")
names(temp)
colnames(df) <- colnames(temp)
rm(temp)

head(df[,1:10])
df$sasdate


df.xts <- xts(df[, -1], as.Date(df$sasdate, "%m/%d/%Y"))
colnames(df.xts)
df.xts <- na.locf(df.xts)

df.xts$S.P.500

plot(df.xts$S.P.500, log(df.xts))


getWeights <- function(d, N)
{
  transToRArray <- function (k)
  {
    k + 1
  }
  w <- rep(NA, N)
  #kth lag
  k = 0
  w[transToRArray(k)] = 1
  for(k in 1:(N-1))
  {
    w[transToRArray(k)] = -1 * w[transToRArray(k-1)] * ( d - k + 1)/k
  }
  w
}

getWeights(0.01, 3)
getWeights(2.5, 3)

#-------------------------------------------------------#
#-----  Replicate 5.1
#-------------------------------------------------------#
N = 6
ds <- seq(from = 0.0, to = 1.0, by = 0.25)
temp <- as.data.frame( sapply(ds, function(d) getWeights(d,N)))
colnames(temp) <- as.character(ds)
temp[["lags"]] <- 0:(N-1)
temp.melt <- reshape2::melt(temp, id.vars = "lags")
head(temp.melt)

library(ggplot2)
ggplot2::ggplot(data = temp.melt, aes(x = lags, y = value, group = variable, colour = variable)) + 
  geom_point( aes(shape = variable)) + 
  geom_line( aes(linetype = variable))



#-------------------------------------------------------#
#-----  Replicate 5.2
#-------------------------------------------------------#
ds <- seq(from = 1, to = 2.0, by = 0.25)
temp <- as.data.frame( sapply(ds, function(d) getWeights(d,N)))
colnames(temp) <- as.character(ds)
temp[["lags"]] <- 0:(N-1)
temp.melt <- reshape2::melt(temp, id.vars = "lags")
head(temp.melt)

library(ggplot2)
ggplot2::ggplot(data = temp.melt, aes(x = lags, y = value, group = variable, colour = variable)) + 
  geom_point( aes(shape = variable)) + 
  geom_line( aes(linetype = variable))


#-------------------------------------------------------#
#-----  FracDiff Expanding Window
#For thresh = 1, nothing skipped
#    d: Any positive fractional
#-------------------------------------------------------#
series = df.xts['1998/2015', "S.P.500"]
#d = 0.4; 
#thresh = 1
#thresh = 0.01
fracdiff.mine <- function(series, d, thresh)
{
  w_       <- getWeights(d, N = length(series))
  rev.w_   <- rev(w_)
  cums     <- cumsum(abs(rev.w_))/sum(abs(rev.w_))
  #This is length of weights vector out of infinity which will be used
  #len.wts  <- length(cums) - max(which(cums < thresh))
  #---------------------------------------------
  #Alternatively we can set all weights below this max index equal to 0
  rev.w_[cums < thresh] <- 0
  res <- sapply(1:length(series), function(i)
                            {
                              sum(series[1:i] * rev.w_[(length(rev.w_) - i + 1):length(rev.w_)])
                            })
  xts(res, index(series))
}

#Book d = 0.4, thresh  = 1 ==> Mine d = 0.4, thresh = 0
first     <- fracdiff.mine(series, d = 0.4, thresh = 1)
first.res <- merge.xts(series, first)
plot(first.res)

#Book d = 0.4, thresh  = 1e-2 ==> Mine d = 0.4, thresh = 0
second <- fracdiff.mine(series, d = 0.4, thresh = 0.001)
second.res <- merge.xts(series, second)
plot(second.res)
tseries::adf.test(second.res$second)	


mytest <- fracdiff.mine(series, d = 1.10, thresh = 0.2)
mytest.res <- merge.xts(series, mytest)
plot(mytest.res)
tseries::adf.test(series)	
tseries::adf.test(mytest.res$mytest)	

#=====================================================================
#-- New Fixed-Width frac diff method
#=====================================================================

getWeights_FFD <- function(d, thresh)
{
  transToRArray <- function (k)
  {
    k + 1
  }
  w <- rep(NA, N)
  #kth lag
  k = 0
  w = 1
  while(1)
  {
    k = k + 1
    wt.temp <- -1 * w[length(w)] * ( d - k + 1)/k
    if(abs(wt.temp) < thresh)
    {
      break
    }
    w = c(w, wt.temp )
  }
  w
}

fracdiff_FFD.mine <- function(series, d, thresh)
{
  w_       <- getWeights_FFD(d, thresh)
  rev.w_   <- rev(w_)
  if(length(rev.w_) < length(series))
  {
    rev.w_ <- c(rep(0, length(series) - length(rev.w_)), rev.w_)
  }
  res <- sapply(1:length(series), function(i)
                                  {
                                    sum(series[1:i] * rev.w_[(length(rev.w_) - i + 1):length(rev.w_)])
                                  })
  xts(res, index(series))
}





mytest     <- fracdiff_FFD.mine(series, d = 0.8, thresh = 0.03)
mytest.res <- merge.xts(series, mytest)
plot(mytest.res)
tseries::adf.test(mytest.res$mytest)	

#Adding a comment



