library(fBasics)
library(xts)
library(PerformanceAnalytics)

getABSMaxDDN <- function(mu, sd, horizon)
{
  #break unit-horizon into n.periods
  n.periods <- 36500
  del.t     <-  1/n.periods
  
  tot.breaks <- horizon * n.periods
  #Dummy but necessary
  #mu*dt + sd*dW
  #dW  N(0,1) * horizon/n.periods
  s.date <- Sys.Date()
  times <- c(s.date, s.date +  c(1:(tot.breaks)))
  
  brownian.shocks <- rnorm(tot.breaks) * sd * sqrt(del.t)
  increments      <- rep(mu * del.t, tot.breaks) + brownian.shocks
  p0 <- 1
  p <- xts(c(1, 1 + cumsum(increments)),times)
  #if(any(p <= 0))
  #{
  #  p[min(which(p <= 0)):length(p)] <- 0  
  #}
  
  #maxddn
  abs.ddns <- as.numeric(cummax(p))  - as.numeric(p[,1])
  #abs.ddns <- abs(log(as.numeric(p[,1])/as.numeric(cummax(p))))
  #plot(abs.ddns, type = "l")
  tail(cummax(abs.ddns), 1)
  #PerformanceAnalytics::Return.calculate(p, "discrete")
}


#Monte carlo
N   = 100000; mu  = 1; sig = 1; horizon = 1
library(parallel)
cl <- makeCluster(4)
clusterExport(cl, c("getABSMaxDDN", "mu", "sig", "horizon"))
ddns <- parallel::parSapply(cl, 1:N, function(i) {
                      library(xts)
                      getABSMaxDDN(mu = mu, sd = sig, horizon = horizon)
        })

stopCluster(cl)

median(ddns)
c(mean(ddns) - qnorm(0.025, lower.tail = FALSE) * sd(ddns)/sqrt(N),  mean(ddns), mean(ddns) + qnorm(0.025, lower.tail = FALSE) * sd(ddns)/sqrt(N))
#1] 0.9172292 0.9194115 0.9215937
hist(ddns)
fBasics::maxddStats(mean = mu, sd = sig, horizon = horizon)
#[1] 0.926318

# %Drawdown:Correct
1- 1/exp(fBasics::maxddStats(mean =  mu - 0.5 * sig * sig, sd = sig , horizon = horizon))
#0.6568503

#Without Drift: Incorrect
1- 1/exp(fBasics::maxddStats(mean =  mu, sd = sig , horizon = horizon))
#0.6039909



