library(manipulate)
library(ggplot2)
lppl <- function(tc, A = 30, B = -1, beta = 0.5, C = 0.1, omega= 20, phi = pi)
{
  t <- seq(from = 0, to = tc, by = 0.001)
  ln_pt = A + (B * (tc-t)^beta) * ( 1 + C * cos(omega * log((tc - t)^beta) + phi))
  dat <- data.frame(time.instant = t, ln_pt)
  ggplot(data = dat, aes(x = time.instant, y = ln_pt) ) + geom_line()
}

manipulate(lppl(tc = tc, A),  
           tc=slider(100,500, initial = 150, step = 50),
           #This is the price at critical time
           A = slider(1, max = 100, initial = 30, step = 1),
           #B < 0: Must be
           B = slider(min = -10, max = 1, initial = -1, step = 1),
           #C : |C| < 1 This param controls magnitude of oscillations
           C  = slider(min = -0.99, max = 0.99, initial = 0.1, step = 0.01)
           )



