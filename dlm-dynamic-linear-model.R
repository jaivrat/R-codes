library(reshape2)
library(ggplot2)
library(dlm)


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
y      <- sum(df[1, ] * c(theta1[1], theta2[1])) + rnorm(1, mean = 0, sd = 0.01)

for(i in 2:1000)
{
  theta1 <- c(theta1, theta1[i-1] + rnorm(1, mean = 0, sd = 0.01))
  theta2 <- c(theta2, theta2[i-1] + rnorm(1, mean = 0, sd = 0.02))
  y      <- c(y, sum(df[i, ] * c(theta1[i], theta2[i])) + rnorm(1, mean = 0, sd = 0.01))
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
res.filt$m

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



library(ggpubr)
ggqqplot(res.fwd.filt$f - df$y[-(1:500)])
plot(res.fwd.filt$f, df$y[-(1:500)])
