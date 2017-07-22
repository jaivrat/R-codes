
N <- 1000


#Cholesky decomposition
e1 <- rnorm(N)
e2 <- rnorm(N)
rho <- 0.3
x1 <- e1
x2 <- rho * x1 + sqrt(1-rho*rho) * e2


dat = data.frame(x = x1, y = x2, type=as.character(rho))

head(dat)

set.seed(1)
e1 <- rnorm(N)
e2 <- rnorm(N)
rho <- 0.8
x1 <- e1
x2 <- rho * x1 + sqrt(1-rho*rho) * e2

dat = rbind(data.frame(x = x1, y = x2, type=as.character(rho)), 
                       dat)

library(ggplot2)

# Set color by cond
ggplot(dat, aes(x=x, y=y, color=type)) + geom_point(shape=1) + 
  geom_smooth(method=lm,   # Add linear regression lines
              se=FALSE)    # Don't add shaded confidence region


#====================== SVD and root matrix way ==============#
covMat <- rbind(c(1, 0.9),
                c(0.9, 1))

svdRes <- svd(covMat)
U <- svdRes$u
D <- diag(svdRes$d)
V <- svdRes$v

#We should get back covMat here
U %*% D %*% t(V)

#Since covMat symmetric => 
eigen(covMat)
A = U %*% D^0.5


getTransformMatrix <- function(covMat)
{
  svdRes <- svd(covMat)
  
  U <- svdRes$u
  D <- diag(svdRes$d)
  V <- svdRes$v
  A = U %*% (D^0.5)
  A
}

rho = 0.8
covMat = rbind(c(1, rho),
               c(rho, 1))

A = getTransformMatrix(covMat)



#First, Cholesky way
#set.seed(1)
e1 <- rnorm(N)
e2 <- rnorm(N)
rho <- 0.8
x1 <- e1
x2 <- rho * x1 + sqrt(1-rho*rho) * e2


#Second spectral transform way
covMat = rbind(c(1, rho),
               c(rho, 1))

A = getTransformMatrix(covMat)
head(data.frame(x1, x2))
head(e1)
head(e2)
E <- rbind(e1, e2)
result <- A %*% E

dat_s <- data.frame(x = result[1,], y = result[2,], "type" = "spectral")
dat_c <- data.frame(x = x1 , y = x2, "type" = "cholesky")
# Set color by cond
ggplot(rbind(dat_s, dat_c), aes(x=x, y=y, color=type)) + geom_point(shape=1) + 
  geom_smooth(method=lm,   # Add linear regression lines
              se=TRUE)    # Don't add shaded confidence region



#------------  Two correlated sets of normals ----------#

rho1 <- 0.8
rho2 = -0.4
#1. Generate independent Normals
e <- matrix(rnorm(2*N), nrow = 2)

covMat1 <- rbind(c(1, rho1),
                c(rho1, 1))
covMat2 <- rbind(c(1, rho2),
                 c(rho2, 1))

result1 <- getTransformMatrix(covMat1) %*% e
result2 <- getTransformMatrix(covMat2) %*% e

dat_1 <- data.frame(x = result1[1,], y = result1[2,], "type" = as.character(rho1))
dat_2 <- data.frame(x = result2[1,], y = result2[2,], "type" = as.character(rho2))

ggplot(rbind(dat_1, dat_2), aes(x=x, y=y, color=type)) + geom_point(shape=1) + 
  geom_smooth(method=lm,   # Add linear regression lines
              se=FALSE)    # Don't add shaded confidence region


