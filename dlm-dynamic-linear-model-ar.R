library(reshape2)
library(ggplot2)
library(dlm)
library(ggpubr)
library(gridExtra)

#-------------------------------------------------------------------------------
#Generate some data first - toy example
#Kalman filter/dynamic linear model in R
#-------------------------------------------------------------------------------
rm(list=ls())
set.seed(100)
df <- data.frame(x1 = sample(seq(from = -1.0, by = 0.01, length.out = 1000), 
                             size = 1000, replace = FALSE),
                 x2 = sample(seq(from = -3.0, by = 0.01, length.out = 1000), 
                             size = 1000, replace = FALSE))

theta1 <- 0.4
theta2 <- 1.2

er1 <- 0.01; phi1 = 0.5
er2 <- 0.05; phi2 = -0.25
err    <- phi1 * er1 + phi2 * er2 + rnorm(1, mean = 0, sd = 0.1)
y      <- sum(df[1, ] * c(theta1[1], theta2[1])) + err
er2 <- er1
er1 <- err

for(i in 2:1000)
{
  theta1 <- c(theta1, theta1[i-1] + rnorm(1, mean = 0, sd = 0.01))
  theta2 <- c(theta2, theta2[i-1] + rnorm(1, mean = 0, sd = 0.02))
  #theta1 <- c(theta1, theta1[i-1])
  #theta2 <- c(theta2, theta2[i-1])
  err    <- phi1 * er1 + phi2 * er2 + rnorm(1, mean = 0, sd = 0.1)
  y      <- c(y, sum(df[i, ] * c(theta1[i], theta2[i])) + err)
  #propagate
  er2 <- er1
  er1 <- err
}


df[["y"]] <- y
df[["tim"]] <- 1:length(y)



#-------------------------------------------------------------------------------
#Learn MLE from first 500
#-------------------------------------------------------------------------------

plot(df[, c("x1", "x2")])

buildFn <- function(x)
{
  ret.mod <- dlm::dlmModReg(X = as.matrix(df[1:500, c("x1", "x2")]),
                            dV = exp(x[3]), 
                            dW = exp(x[1:2]), 
                            addInt = FALSE)
  ret.mod
}

fit <- dlm::dlmMLE( y     = df$y[1:500], 
                    build = buildFn, 
                    parm  = log(c(rep(1e2, 2), 1e2)),
                    lower = log(rep(1e-5, 3)), 
                    hessian=TRUE)


mod.res <- buildFn(fit$par)
avarLog <- solve(fit$hessian)
avar <- diag(exp(fit$par)) %*% avarLog %*% diag(exp(fit$par)) # Delta method
V(mod.res)
W(mod.res)


res.filt <- dlm::dlmFilter(y = df$y[1:500], mod = mod.res)
tail(res.filt$m)

results.df <- data.frame(res.filt$m[-1, ] , 
                         theta1 = theta1[1:500], 
                         theta2 = theta2[1:500],
                         tim = df$tim[1:500])
head(results.df)

library(reshape2)
results.df.melt <- reshape2::melt(results.df, id.vars = "tim")     
ggplot(data = results.df.melt, aes(x = tim, y = value, group = variable, colour = variable)) + 
  geom_line() 


#-------------------------------------------------------------------------------
#Future prediction though kalman filter
#-------------------------------------------------------------------------------
ret.mod.fwd <- dlm::dlmModReg(X = as.matrix(df[-(1:500), c("x1", "x2")]),
                              dV = V(mod.res), 
                              dW = diag(W(mod.res)), 
                              m0 = as.numeric(tail(res.filt$m, 1)),
                              C0 = with(res.filt, dlmSvd2var(U.C[[501]], D.C[501,])),
                              addInt = FALSE)

res.fwd.filt   <- dlm::dlmFilter(y = df$y[-(1:500)], mod = ret.mod.fwd)

results.fwd.df <- data.frame(res.fwd.filt$m[-1, ] , 
                             theta1 = theta1[-(1:500)], 
                             theta2 = theta2[-(1:500)],
                             tim = df$tim[-(1:500)])
head(results.fwd.df)


results.fwd.df.melt <- reshape2::melt(results.fwd.df, id.vars = "tim")     
ggplot(data = results.fwd.df.melt, aes(x = tim, y = value, group = variable, colour = variable)) + 
  geom_line() + 
  ggtitle("States, ie thetas")


#plot y's and predictions
results.y.df  <- data.frame(tim = df$tim[-(1:500)], y = df$y[-(1:500)], pred = res.fwd.filt$f)
results.y.df.melt <- reshape2::melt(results.y.df, id.vars = "tim")     
ggplot(data = results.y.df.melt, aes(x = tim, y = value, group = variable, colour = variable)) + 
  geom_line() + 
  ggtitle("ys")

#Normality plot
ggqqplot(res.fwd.filt$f - df$y[-(1:500)])
plot(res.fwd.filt$f, df$y[-(1:500)])

acf(res.fwd.filt$f - df$y[-(1:500)])
head(residuals(res.fwd.filt, type = "raw")$res)
head(res.fwd.filt$f - df$y[-(1:500)])

head(scale(residuals(res.fwd.filt, type = "raw")$res))
head(scale(res.fwd.filt$f - df$y[-(1:500)]))

head(residuals(res.fwd.filt, type = "standardized")$res)
head((res.fwd.filt$f - df$y[-(1:500)])/sd(res.fwd.filt$f - df$y[-(1:500)]))

sd(res.fwd.filt$f - df$y[-(1:500)])
#residuals(res.fwd.filt)$sd
with(res.fwd.filt, dlmSvd2var(U.R[[500]], D.R[500,]))


#what is sd term , and how it is related to the filter above
tail(residuals(res.fwd.filt)$sd, 1)^2
#Formula Directly from algo Q_t = F_t %*% R %*% F_t' + V
tail(res.fwd.filt$mod$X, 1) %*% with(res.fwd.filt, dlmSvd2var(U.R[[500]], D.R[500,])) %*% t(tail(res.fwd.filt$mod$X, 1)) + res.fwd.filt$mod$V




#---------------------------------------------------------------------------------
#AS WE SEE ABOVE ACF of plor shows autocorrelations
#---------------------------------------------------------------------------------
acf(res.fwd.filt$f - df$y[-(1:500)])

#We expect it to be AR(2) as we have generated data that way
#Build new function with AR components
#parm  = log(c(rep(1e2, 2), 1e2)),
#lower = log(rep(1e-5, 3)), 
ar.order = 3
buildFn.ar <- function(x)
{
  ret.mod <- dlm::dlmModReg(X = as.matrix(df[1:500, c("x1", "x2")]),
                            dV = exp(x[3]), 
                            dW = exp(x[1:2]), 
                            addInt = FALSE)
  #ar: ar coefficients
  ret.arma <- dlm::dlmModARMA(ar = x[4:(4 + ar.order-1)], sigma2 =  exp(x[4 + ar.order]))
  return(ret.mod + ret.arma )
}


#mod.ar <- buildFn.ar(log(c(rep(1e2, 2), 1e2)))
fit.ar <- dlm::dlmMLE( y     = df$y[1:500], 
                    build = buildFn.ar, 
                    parm  = c( log(c(rep(1, 2), 1e2))
                               ,rep(1, ar.order) #ar components
                               ,log(1) #sigma for ar
                               ),
                    lower = c( log(rep(1e-5, 3))
                                #ar compoenents
                                ,rep(-Inf, ar.order)
                                ,log(1e-5)
                                ), 
                    hessian=TRUE)
mod.ar.res <- buildFn.ar(fit.ar$par)
avarLog <- solve(fit.ar$hessian)
avar <- diag(exp(fit.ar$par)) %*% avarLog %*% diag(exp(fit.ar$par)) # Delta method
V(mod.ar.res)
W(mod.ar.res)


#Filter to determine params
res.ar.filt <- dlm::dlmFilter(y = df$y[1:500], mod = mod.ar.res)
res.ar.filt$m


#-------------------------------------------------------------------------------
#Future prediction though kalman filter
#-------------------------------------------------------------------------------
buildFn.ar.fwd <- function(x)
{
  ret.mod <- dlm::dlmModReg(X = as.matrix(df[-(1:500), c("x1", "x2")]),
                            dV = exp(x[3]), 
                            dW = exp(x[1:2]), 
                            addInt = FALSE)
  #ar: ar coefficients
  ret.arma <- dlm::dlmModARMA(ar = x[4:(4 + ar.order-1)], sigma2 =  exp(x[4 + ar.order]))
  return(ret.mod + ret.arma )
}

ret.mod.ar.fwd <- buildFn.ar.fwd(fit.ar$par)
#replace m and C0 with filder
ret.mod.ar.fwd$m0
ret.mod.ar.fwd$m0 <- tail(res.ar.filt$m, 1)
ret.mod.ar.fwd$C0 <- with(res.ar.filt, dlmSvd2var(U.C[[length(res.ar.filt$U.C)]], D.C[length(res.ar.filt$U.C),]))
#We can get phi1 and phi2 from this
ret.mod.ar.fwd$GG


res.fwd.ar.filt   <- dlm::dlmFilter(y = df$y[-(1:500)], mod = ret.mod.ar.fwd)

results.fwd.ar.df <- data.frame(res.fwd.ar.filt$m[-1, 1:2], 
                                theta1 = theta1[-(1:500)], 
                                theta2 = theta2[-(1:500)],
                                tim = df$tim[-(1:500)])
head(results.fwd.ar.df)
results.fwd.df.ar.melt <- reshape2::melt(results.fwd.ar.df, id.vars = "tim")     
ggplot(data = results.fwd.df.ar.melt, aes(x = tim, y = value, group = variable, colour = variable)) + 
  geom_line() + 
  ggtitle("States, ie thetas")

#plot y's and predictions
results.y.ar.df  <- data.frame(tim = df$tim[-(1:500)], y = df$y[-(1:500)], pred = res.fwd.ar.filt$f)
results.y.ar.df.melt <- reshape2::melt(results.y.ar.df, id.vars = "tim")     
ggplot(data = results.y.ar.df.melt, aes(x = tim, y = value, group = variable, colour = variable)) + 
  geom_line() + 
  ggtitle("ys")

#Normality plot
ggqqplot(res.fwd.ar.filt$f - df$y[-(1:500)])
plot(res.fwd.ar.filt$f, df$y[-(1:500)])

head(residuals(res.fwd.ar.filt, type = "raw")$res)
head(res.fwd.ar.filt$f - df$y[-(1:500)])

head(scale(residuals(res.fwd.ar.filt, type = "raw")$res))
head(scale(res.fwd.ar.filt$f - df$y[-(1:500)]))

head(residuals(res.fwd.ar.filt, type = "standardized")$res)
head((res.fwd.ar.filt$f - df$y[-(1:500)])/sd(res.fwd.ar.filt$f - df$y[-(1:500)]))

sd(res.fwd.ar.filt$f - df$y[-(1:500)])
#residuals(res.fwd.filt)$sd
with(res.fwd.ar.filt, dlmSvd2var(U.R[[500]], D.R[500,]))


#what is sd term , and how it is related to the filter above
tail(residuals(res.fwd.ar.filt)$sd, 1)^2
#Formula Directly from algo Q_t = F_t %*% R %*% F_t' + V
#tail(res.fwd.ar.filt$mod$X, 1) %*% with(res.fwd.ar.filt, dlmSvd2var(U.R[[500]], D.R[500,])) %*% t(tail(res.fwd.ar.filt$mod$X, 1)) + res.fwd.ar.filt$mod$V


shapiro.test(res.fwd.ar.filt$f - df$y[-(1:500)])
#TEST does not still show normality :(  Oh! My god - what to do!! - need a math professor now. I feel I should weep -Jai
#Call.. lord Ram/Hanuman-ji

#Three necessary plots
ggqqplot(res.fwd.ar.filt$f - df$y[-(1:500)])
plot(res.fwd.ar.filt$f, df$y[-(1:500)])
acf(residuals(res.fwd.ar.filt, type = "raw")$res)

