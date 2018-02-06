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
    K * Q1 /r
  }
  v2 <- function(x,y)
  {
    r <- sqrt((D-x)*(D-x) + y*y)
    K * Q2 /r
  }
  
  v1(x,y) + v2(x,y)
}


E <- function(x,y, Q1, Q2)
{
   eps <- 1e-6
   E_x = -(V(x + eps, y, Q1, Q2) - V(x - eps, y, Q1, Q2))/(2 * eps)
   E_y = -(V(x, y + eps, Q1, Q2) - V(x, y - eps, Q1, Q2))/(2 * eps)
   list(E_x=E_x, E_y=E_y)
}

#Directions make sense
E(4,0, Q1, Q2)
E(9,0, Q1, Q2)
E(-4,0, Q1, Q2)
E(0,4, Q1, Q2)

#at what theta, the angle is 60deg
theta1 = 0
theta2 = pi/2

bisection <- function(a,b, angle)
{
  obj   = tan(angle)
  
  if( (tan(b) - obj) * (tan(a) - obj) > 0)
    stop("Wrong choice of intervals")
  
  mn = min(a,b)
  mx = max(a,b)
  a = mn; b = mx
  #initial case
  theta = (a + b)/2
  if(abs(tan(theta) - obj) < 1e-10)
  {
    return(theta)
  }
  
  while(abs(tan(theta) - obj) > 1e-10)
  {
    print(sprintf("a=%s b=%s, theta = %s\n", a*180/pi, b*180/pi, theta * 180/pi))
    if( (tan(theta) - obj) * (tan(a) - obj) < 0 ){
      b = theta
    } else if ( (tan(theta) - obj) * (tan(b) - obj) < 0 ) {
      a = theta
    }
    theta = (a + b)/2
  }
  return(theta)
}

#At this point angle is 60
theta = bisection(0, pi/2, angle = pi/3)
print(sprintf("Angle = %s", theta * 180/pi))
x0 = R * cos(theta)
y0 = R * sin(theta)
eField = E(x0,y0, Q1, Q2)
print(eField)
print(sprintf("Angle of field = %s degrees" , atan(eField$E_y/eField$E_x) * 180/pi))


#from this we find the curve











