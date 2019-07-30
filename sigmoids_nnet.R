sigmoid_and <- function(x)
{
  theta.t.x <- sum(c(1,x) * c(-30, 20, 20))
  1/(1+exp(-theta.t.x))
}

round(sigmoid_and(c(0,0)), 2)
round(sigmoid_and(c(0,1)), 2)
round(sigmoid_and(c(1,0)), 2)
round(sigmoid_and(c(1,1)), 2)

sigmoid_or <- function(x)
{
  theta.t.x <- sum(c(1,x) * c(-10, 20, 20))
  1/(1+exp(-theta.t.x))
}

round(sigmoid_or(c(0,0)), 2)
round(sigmoid_or(c(0,1)), 2)
round(sigmoid_or(c(1,0)), 2)
round(sigmoid_or(c(1,1)), 2)



sigmoid_not <- function(x)
{
  theta.t.x <- sum(c(1,x) * c(10, -20))
  1/(1+exp(-theta.t.x))
}

round(sigmoid_not(c(0)), 2)
round(sigmoid_not(c(1)), 2)


#-----------------------------------------
#XNOR
#-----------------------------------------
nnet_xnor <- function(x)
{
  a1 <- sigmoid_and(x)
  a2 <- 1/(1+exp(-sum(c(10,-20,-20) * c(1,x))))
  sigmoid_or(x = c(a1,a2))
}

round(nnet_xnor(c(0,0)), 2)
round(nnet_xnor(c(0,1)), 2)
round(nnet_xnor(c(1,0)), 2)
round(nnet_xnor(c(1,1)), 2)


