rm(list=ls())
library(ggplot2)

D = 5
R = 0.5
epsilon <- 8.85e-12
#pi
#K = 1/(4*pi*epsilon)
K = 1
Q <- 1
Q1 <- 4*Q
Q2 <- -Q
#Q2 <- -4*Q

V_grad <- function(x,y, Q1, Q2)
{
  v1_grad <- function(x,y)
  {
    list(dvx = - K * Q1 * x * (x^2 + y^2)^(-3/2), dvy = - K * Q1 * y * (x^2 + y^2)^(-3/2))
  }
  
  v2_grad <- function(x,y)
  {
    list(dvx = K * Q2 * (D - x) * ((D-x)^2 + y^2)^(-3/2),  dvy = -K * Q2 * y *  ((D-x)^2 + y^2)^(-3/2))
  }
  
  v1_grad_res = v1_grad(x,y)
  v2_grad_res = v2_grad(x,y)
  
  list(V_grad_x = v1_grad_res$dvx + v2_grad_res$dvx,  V_grad_y = v1_grad_res$dvy + v2_grad_res$dvy)
}
  
E <- function(x,y, Q1, Q2)
{
   grads = V_grad(x,y, Q1, Q2)
   list(E_x=-grads$V_grad_x, E_y=-grads$V_grad_y)
}


#Directions make sense
E(4,0, Q1, Q2); #E.old(4,0, Q1, Q2)
E(9,0, Q1, Q2); #E.old(9,0, Q1, Q2)
E(-4,0, Q1, Q2); #E.old(-4,0, Q1, Q2)
E(0,4, Q1, Q2)

#at what theta, the angle is 60deg
theta1 = 0
theta2 = pi/2

bisection <- function(a,b, angle)
{
  obj   = tan(angle)
  
  if( (tan(b) - obj) * (tan(a) - obj) > 0)
    stop("Wrong choice of intervals")
  
  if(tan(a) == obj)
  {
    return(a)
  }
  
  if(tan(b) == obj)
  {
    return(b)
  }
  
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
    print(sprintf("bisection: a=%s b=%s, theta = %s\n", a*180/pi, b*180/pi, theta * 180/pi))
    if( (tan(theta) - obj) * (tan(a) - obj) < 0 ){
      b = theta
    } else if ( (tan(theta) - obj) * (tan(b) - obj) < 0 ) {
      a = theta
    }
    theta = (a + b)/2
  }
  return(theta)
}


#angle = pi/3
getFieldPoints <- function(angle)
{
  #At this point angle is 60
  theta = bisection(0, pi/2, angle)
  print(sprintf("getFieldPoints: Angle = %s", theta * 180/pi))
  x0 = R * cos(theta)
  y0 = R * sin(theta)
  eField = E(x0,y0, Q1, Q2)
  print(eField)
  print(sprintf("getFieldPoints: Angle of field = %s degrees" , atan(eField$E_y/eField$E_x) * 180/pi))
  
  #from this we find the curve
  getNewCord <- function(x,y)
  {
    ds <- 0.001
    eField <- E(x,y, Q1, Q2)
    theta <- abs(atan(eField$E_y/eField$E_x))
    x_new <- x + ds * cos(theta) * sign(eField$E_x)
    y_new <- y + ds * sin(theta) * sign(eField$E_y)
    return(c(x_new, y_new))
  }
  
  CheckEntry <- function(x,y)
  {
    ((x - D)^2 + (y -0)^2 < R*R) || (x >= D + 30*R)
  }
  
  getPoints <- function(x0,y0)
  {
    entered = FALSE
    x = x0; y = y0;
    df <- data.frame(x=x0,y=y0)
    while(!entered)
    {
      print(sprintf("getPoints(angle=%s): X=%s, Y=%s\n", round(angle*180/pi,2), x,y))
      pt <- getNewCord(x,y)
      df <- rbind(df,c(pt[1], pt[2]))
      x = pt[1]; y = pt[2]
      entered = CheckEntry(x, y)
    }
    df
  }
  
  #debug(getPoints)
  ptdf <- getPoints(x0,y0)
  ptdf[["TYPE"]] <- sprintf("%s", round(angle * 180/pi,2))
  
  ptdf
}


dataDF <- rbind(getFieldPoints(0),
                getFieldPoints(pi/6),
                getFieldPoints(pi/4),
                getFieldPoints(pi/3))

# Extend the regression lines beyond the domain of the data
ggplot(dataDF, aes(x=x, y=y, color=TYPE)) + geom_point(shape=1) +
  scale_colour_hue(l=50)  # Use a slightly darker palette than normal


#add another circle
getCircle <- function(radius, center)
{
  xs <- seq(from = center[1] - radius, to = center[1] + radius, length.out = 1000)
  ys <- sqrt( radius*radius - (xs - center[1])^2)
  ysPos <- ys + center[2]
  ysNeg <- -ys + center[2]
  posDF <- data.frame(x = xs, y = ysPos)
  negDF <- data.frame(x = xs, y = ysNeg)
  return(rbind(posDF, negDF))
}


circdf1 <- getCircle(R, c(0,0))
circdf1[["TYPE"]] <- "sphere1"
circdf2 <- getCircle(R, c(D,0))
circdf2[["TYPE"]] <- "sphere2"

allDF <- rbind(dataDF, circdf1, circdf2)



# Extend the regression lines beyond the domain of the data
ggplot(allDF, aes(x=x, y=y, color=TYPE)) + geom_point(shape='....') +
  scale_colour_hue(l=50)   # Use a slightly darker palette than normal


#Specifically for 60deg : When charge equal
#debug(getFieldPoints)
df60 <- getFieldPoints(pi/4)
lastpoint <- df60[dim(df60)[1],]
fields <- E(lastpoint$x,lastpoint$y, Q1, Q2)
atan(abs(fields$E_y/fields$E_x)) * 180/pi

# Extend the regression lines beyond the domain of the data
ggplot(df60, aes(x=x, y=y, color=TYPE)) + geom_point(shape='....') +
  scale_colour_hue(l=50)   # Use a slightly darker palette than normal
