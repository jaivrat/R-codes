#This is not a serious model - however to go through a cycle of usage on a standard dataset.
rm(list=ls())
library(PerformanceAnalytics)
library(dlm)
library(reshape2)
library(ggplot2)
library(ggpubr)
library(gridExtra)
library(xts)

df <- read.csv(file = "/Users/jvsingh/work/github/rcodes/2019-04.txt", 
               stringsAsFactors = FALSE, check.names = FALSE)
df <- df[-1,]
#head(df)
df$sasdate <- as.Date(df$sasdate, "%m/%d/%Y")
head(df[, 1:10])



sp.price.xts <-  xts(df$`S&P 500`, df$sasdate)
sp.log.ret.xts <- PerformanceAnalytics::Return.calculate(sp.price.xts, method = "log")
sp.sim.ret.xts <- PerformanceAnalytics::Return.calculate(sp.price.xts, method = "discrete")

other.dat.xts  <- xts(df[, 2:3], df$sasdate)

in.ret.xts      <- sp.log.ret.xts
#---function : calculateAdjustedReturns
calculateAdjustedReturns <- function(in.ret.xts,
                                     i.start.period,
                                     j.span.periods)
{
   roll.ret <- PerformanceAnalytics::apply.rolling(in.ret.xts, 
                                                   width = j.span.periods,
                                                   trim  = FALSE, 
                                                   FUN   = sum)
   shited.xts <- xts::lag.xts(roll.ret, 
                              k      = -(i.start.period + j.span.periods), 
                              na.pad = TRUE )
   res.xts    <- xts::merge.xts(in.ret.xts, roll.ret, shited.xts)
   names(res.xts) <- c("in.ret", "roll.ret", "shifted.ret")
   res.xts
}

my.dlmForecast <- function(m.start, C.start, G, fut.FF, W, V, k.steps)
{
  a.list <- list()
  R.list <- list()
  f.list <- list()
  Q.list <- list()
  
  a.list[["0"]] <- matrix(m.start, ncol = 1)
  R.list[["0"]] <- C.start
  for(k in 1:k.steps)
  {
    #Theta
    a.list[[as.character(k)]] <- G %*% a.list[[as.character(k-1)]]
    R.list[[as.character(k)]] <- G %*% R.list[[as.character(k-1)]] %*% t(G) + W
    
    #Y
    f.list[[as.character(k)]] <- fut.FF[k, ] %*% a.list[[as.character(k)]]
    Q.list[[as.character(k)]] <- fut.FF[k, ] %*% R.list[[as.character(k)]] %*% t(fut.FF[k, ]) + V
  }
  
  list(a = a.list[[length(a.list)]], 
       R = R.list[[length(R.list)]], 
       f = as.numeric(f.list[[length(f.list)]]), 
       Q = as.numeric(Q.list[[length(Q.list)]]))
}

i.start.period =0
j.span.periods =2
adj.xts <- calculateAdjustedReturns(in.ret.xts = sp.log.ret.xts,
                                    i.start.period = i.start.period,
                                    j.span.periods = j.span.periods)
#################################################################################
#--lets apply
dat.xts <- merge(other.dat.xts, adj.xts)
dat.xts <- dat.xts[, c("RPI", "W875RX1", "shifted.ret")]
dat.xts <- dat.xts[-1, ] #First is NA
#tail(dat.xts, 10)
#tail(other.dat.xts[,1:2], 10)
#tail(sp.price.xts)
#dat.xts <- dat.xts[complete.cases(dat.xts),]
dat.mle <- dat.xts["1900/1998"]

buildFn <- function(x)
{
  ret.mod <- dlm::dlmModReg(X = as.matrix(dat.mle[,1:2]),
                            dV = exp(x[3]), 
                            dW = exp(x[1:2]), 
                            addInt = FALSE)
  ret.mod
}

fit <- dlm::dlmMLE( y     = dat.mle$shifted.ret, 
                    build = buildFn, 
                    parm  = log(c(rep(1e2, 2), 1e2)),
                    lower = log(rep(1e-5, 3)), 
                    hessian=TRUE)
mod.res <- buildFn(fit$par)
avarLog <- solve(fit$hessian)
avar <- diag(exp(fit$par)) %*% avarLog %*% diag(exp(fit$par)) # Delta method
V(mod.res)
W(mod.res)


res.filt <- dlm::dlmFilter(y = dat.mle$shifted.ret, mod = mod.res)
res.filt$m
plot.xts(xts(res.filt$m[-1,], index(dat.mle)))


#-------------------------------------------------------------------------------
#Future prediction though kalman filter
#-------------------------------------------------------------------------------
dat.test  <- dat.xts["1999/3000"]
curr.date <- as.Date(index(dat.test))[1]
getPred.RealScenario <- function(curr.date)
{
  dat.filt.xts <- dat.xts[paste("1900", "/" , as.character(curr.date), sep= "")]
  i.plus.j     <- i.start.period  + j.span.periods
  
  #Complete realised y's in this
  dat.filt.realised.xts <- dat.filt.xts[1:(dim(dat.filt.xts)[1]-i.plus.j)]
  #Rest of the y's are blank
  dat.filt.unrealised.xts <- dat.filt.xts[(dim(dat.filt.xts)[1]-i.plus.j + 1):(dim(dat.filt.xts)[1])]
  
  #We run filter up to above to determine the latest m0 and C0
  ret.mod.fwd <- dlm::dlmModReg(X = as.matrix(dat.filt.realised.xts[,1:2]),
                                dV = V(mod.res), 
                                dW = diag(W(mod.res)),
                                m0 = rep(0, 2),
                                C0 = diag(1e7, 2),
                                addInt = FALSE)
  res.fwd.filt   <- dlm::dlmFilter(y = dat.filt.realised.xts$shifted.ret, mod = ret.mod.fwd)
  
  #Using last m0 and C0 simulate for forward unlrealised, in which the last one will 
  #give the required y metric of ineterest 
  m0 = tail(res.fwd.filt$m,1)
  C0 = with(res.fwd.filt, dlmSvd2var(U.C[[length(res.fwd.filt$U.C)]], D.C[length(res.fwd.filt$U.C),]))
  
  res <- my.dlmForecast(m.start = m0, C.start =  C0, G = diag(2), 
                        fut.FF = dat.filt.unrealised.xts[,1:2], 
                        W = ret.mod.fwd$W, 
                        V = ret.mod.fwd$V, 
                        k.steps = i.plus.j)
  
  c(res$f, res$Q)
}

#These results are such that for each date we have i+j missing returns 
#because we are using future data at each stamp.
#at (1) => T, (2) =>T-1, ..(i+j) => T-(i+j-1)
final.results <- sapply(as.Date(index(dat.test)), getPred.RealScenario)

#Lets us check how it went
result.xts <- xts(data.frame(t(final.results)), as.Date(index(dat.test)))
names(result.xts) <- c("f", "Q")

plot.xts(result.xts$f)
plot.xts(result.xts$Q)

tmp <- merge.xts(dat.xts$shifted.ret, result.xts$f)
plot.xts(tmp)
plot.xts(tmp[,1] - tmp[,2])

#Other Check - if we directly run filter on all the data and see last values
ret.mod.all <- dlm::dlmModReg(X = as.matrix(dat.xts[,1:2]),
                              dV = V(mod.res), 
                              dW = diag(W(mod.res)),
                              m0 = rep(0, 2),
                              C0 = diag(1e7, 2),
                              addInt = FALSE)
res.filt.all   <- dlm::dlmFilter(y = dat.xts$shifted.ret, mod = ret.mod.all)


#When i=0, j = 1. Do not expect difference here:
#   Because both cases use same theta to determine the one-step forecast.
#When i=0, j = 2. Expect difference here. 
#   Because result.xts uses two period old theta to predit 2 period later theta
tail(res.filt.all$f)
tail(result.xts$f)

