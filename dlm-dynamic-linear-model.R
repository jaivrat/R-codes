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

#Normality plot
ggqqplot(res.fwd.filt$f - df$y[-(1:500)])
plot(res.fwd.filt$f, df$y[-(1:500)])

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
#Let us check if MLE estimates of V and W matter
#---------------------------------------------------------------------------------
#with(res.fwd.filt, dlmSvd2var(U.R, D.R))
#with(res.fwd.filt, dlmSvd2var(U.R[[500]], D.R[500,]))
#with(res.fwd.filt, dlmSvd2var(U.R[[499]], D.R[499,]))
residuals.sd.true <- residuals(res.fwd.filt)$sd ^2

#model on some choice of dV and dW
V(mod.res)
#[1,] 0.0001028086
W(mod.res)
# > W(mod.res)
#              [,1]       [,2]
# [1,] 0.0001065652 0.00000000
# [2,] 0.0000000000 0.00038172


library(manipulate)

VisualiseChange <- function(pV.pct.inc, dW.pct.inc)
{
      ##-- all goes in function
      #conctsurct a model, based on assumed dV and dW
      dV.new <- V(mod.res) * ( 1 + pV.pct.inc/100)
      dW.new <- W(mod.res) * ( 1 + dW.pct.inc/100)
      #Use filter on first half to come up with m0 and C0 for next fwd data set.
      mod.first.half <- dlm::dlmModReg(X = as.matrix(df[1:500, c("x1", "x2")]),
                                    dV = dV.new, 
                                    dW = diag(dW.new), 
                                    m0 = rep(0, 2),
                                    C0 = diag(2) * 1e7,
                                    addInt = FALSE)
      res.first.half.filt   <- dlm::dlmFilter(y = df$y[1:500], mod = mod.first.half)
      #Create a new model with assumed dV dW and m0 C0
      m0 = tail(res.first.half.filt$m, 1)
      C0 = with(res.first.half.filt, dlmSvd2var(U.C, D.C))
      mod.second.half <- dlm::dlmModReg(X = as.matrix(df[-(1:500), c("x1", "x2")]),
                                       dV = dV.new, 
                                       dW = diag(dW.new), 
                                       m0 = m0,
                                       C0 = C0[[length(C0)]],
                                       addInt = FALSE)
      #apply filter and see results
      res.second.half.filt   <- dlm::dlmFilter(y = df$y[-(1:500)], mod = mod.second.half)
      
      
      final.df <- data.frame("tim" = df$tim[-(1:500)],
                             "actual.y" = df$y[-(1:500)],
                             "pred.y.true" = res.fwd.filt$f,
                             "pred.y.apprx" = res.second.half.filt$f)
      final.df.melt = reshape2::melt(final.df, "id.vars" = "tim")
      p1 <- ggplot(final.df.melt, aes(x = tim, y = value, group = variable, colour = variable)) + geom_point() + geom_line()
      
      
      
      
      final.df <- data.frame("tim" = df$tim[-(1:500)],
                             "pred.Q" = residuals(res.fwd.filt)$sd,
                             "pred.Q.apprx" = residuals(res.second.half.filt)$sd)
      final.df.melt = reshape2::melt(final.df, "id.vars" = "tim")
      p2 <- ggplot(final.df.melt, aes(x = tim, y = value, group = variable, colour = variable)) + geom_point() + geom_line()
      
      
      #with(res.second.half.filt, dlmSvd2var(U.R, D.R))
      
      var.mag <- data.frame(t(sapply(1:length(df$tim[-(1:500)]), function(i){
                                          Ft = res.fwd.filt$mod$X[i, ]
                                          Qt = dlmSvd2var(res.fwd.filt$U.R[[i]], res.fwd.filt$D.R[i,])
                                          c(Ft %*% Qt %*% Ft, res.fwd.filt$mod$V)
                                        })))
      colnames(var.mag) <- c("QTerm", "V")
      var.mag[["tim"]] = df$tim[-(1:500)]
      var.mag.melt = reshape2::melt(var.mag, "id.vars" = "tim")
      p3 <- ggplot(var.mag.melt, aes(x = tim, y = value, group = variable, colour = variable)) + geom_point() + geom_line() + ggtitle("Qpart & V part")
      
      library(gridExtra)
      p.final <- grid.arrange(
        p1,
        p2,
        p3,
        ncol = 1,
        top = "Title of the page"
      )
      
      ##-- all goes in function
      p.final
} 

manipulate(VisualiseChange(pV.pct.inc = pV.pct.inc, dW.pct.inc = dW.pct.inc), pV.pct.inc = slider(0, 100), dW.pct.inc = slider(0, 100))

#Seems some error in estimation of dV and dW does not matter. Kalman filter process handles it
#errors are more weightage coming out of Q matrix rather than measurement errors V





