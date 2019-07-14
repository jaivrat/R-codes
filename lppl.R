library(manipulate)
library(ggplot2)
lppl <- function(tc, A, B, beta, C, omega, phi)
{
  t <- seq(from = 0, to = tc, by = 0.001)
  ln_pt = A + (B * (tc-t)^beta) * ( 1 + C * cos(omega * log((tc - t)^beta) + phi))
  pt    <- exp(ln_pt)
  dat <- data.frame(time.instant = t, pt, ln_pt)
  ggplot(data = dat, aes(x = time.instant, y = ln_pt) ) + 
    geom_line() + 
    geom_hline(yintercept = A, color = "red")
}

manipulate(lppl(tc = tc, A, B, beta, C, omega, phi = n.of.npi * pi),  
           tc=slider(100,500, initial = 150, step = 50),
           #This is the price at critical time
           A = slider(1, max = 100, initial = 30, step = 1, 
                      label = "A: This is the price at critical time"),
           #B < 0: Must be
           B = slider(min = -10, max = 1, initial = -1, step = 1, 
                      label = "B must be less than 0"),
           #C : |C| < 1 This param controls magnitude of oscillations
           C  = slider(min = -0.99, max = 0.99, initial = 0.1, step = 0.01, 
                       label = "C: |C| < 1 This param controls magnitude of oscillations"),
           #beta: 0 < beta< 1 - controls the growth rate of the magnitude and is the most
           #important feature capturing the imminence of a regime witching, as his value 
           #is close to zero.
           beta = slider(min = 0.01, max = 0.999, initial = 0.5, step = 0.01),
           omega = slider(min = 0, max = 9999, step = 10, initial = 20),
           n.of.npi = slider(min = 0, max = 2 , label = "n*pi", step = 1/2, initial = 1)
           )



