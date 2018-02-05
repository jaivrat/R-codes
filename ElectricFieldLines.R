
D = 10
R = 1
epsilon <- 8.85e-12
#pi
#K = 1/(4*pi*epsilon)
K = 1
Q <- 100
Q1 <- 4*Q
Q2 <- -Q


V <- function(x,y, Q1, Q2)
{
  v1 <- function(x,y)
  {
    r <- sqrt(x*x + y*y)
    K * Q1 * (1/r)
  }
  v2 <- function(x,y)
  {
    r <- sqrt((D-x)*(D-x) + y*y)
    K * Q2 * (1/r)
  }
  
  v1(x,y) + v2(x,y)
}


E <- function(x,y, Q1, Q2)
{
   eps <- 1e-8
   E_x = -(V(x + eps, y, Q1, Q2) - V(x - eps, y, Q1, Q2))/(2 * eps)
   E_y = -(V(x, y + eps, Q1, Q2) - V(x, y - eps, Q1, Q2))/(2 * eps)
   list(E_x, E_y)
}

#Directions make sense
E(4,0, Q1, Q2)
E(9,0, Q1, Q2)
E(-4,0, Q1, Q2)
E(0,4, Q1, Q2)

#at what theta, the abgle is 60deg

theta1 = 0
theta2 = pi











