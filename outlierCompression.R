
library(ggplot2)

#------------------- Function to generate correlated normals(Spectral Transform) -------------------#
getTransformMatrix <- function(covMat)
{
  svdRes <- svd(covMat)
  
  U <- svdRes$u
  D <- diag(svdRes$d)
  V <- svdRes$v
  A = U %*% (D^0.5)
  A
}

dimension=2
N      <- 1000
rho    <- 0.8

covMat <- rbind(c(1, rho),
               c(rho, 1))
A      <- getTransformMatrix(covMat)
e1     <- rnorm(N); e2 <- rnorm(N); E <- rbind(e1, e2)
result <- A %*% E

datf <- data.frame(x = result[1,], y = result[2,], "ORIG" = TRUE) 
#Now let us add some nasty data to datf.
rho    <- 0.8
covMat <- rbind(c(1, rho),
                c(rho, 1))
A      <- getTransformMatrix(rbind(c(6, -0.8),
                                   c(-0.8, 6)))
e1     <- rnorm(50); e2 <- rnorm(50); E <- rbind(e1, e2)
result <- A %*% E
datf.o <- data.frame(x = result[1,], y = result[2,], "ORIG" = FALSE) 
datf <-     rbind(datf,
                  datf.o)

ggplot(datf, aes(x=x, y=y)) + geom_point(shape=1) + 
  geom_smooth(method=lm,   # Add linear regression lines
              se=FALSE)    # Don't add shaded confidence region



findoutliers.mahlanobis <- function(in.df, dimension, pctOuliers)
{
  in.df[["OUTL"]] <- FALSE
  
  #let as fix outlier count
  L = round(pctOuliers * dim(in.df)[1], 0)
  
  m <- matrix()
  delta_cost <- Inf
  i <- 0
  eps <- 1e-10
  prevCov <- NULL
  
  cols <- colnames(in.df)
  while(1)
  {
    i = i + 1
    
    center <- matrix(apply(in.df[in.df$OUTL != TRUE, 1:dimension], 2, mean, na.rm = TRUE), ncol = 1)
    covmat <- cov(in.df[in.df$OUTL != TRUE, 1:dimension])
    
    #Record iteration
    if(i==1)
    {
      m <- matrix(c(i, det(covmat)), nrow = 1)
    } else {
      m <- rbind( m,
                  c(i, det(covmat)))
    }
    
    #mahlanobis distance based on centers and cov
    mh.dist  <- mahalanobis(x = in.df[, 1:dimension], center = center, cov = covmat)
    in.df[["mh.dist"]] <- mh.dist
    in.df <- in.df[order(-mh.dist),]
    
    #discard top L
    in.df[["OUTL"]] <- FALSE
    in.df[1:L, "OUTL"] <- TRUE
    
    newcovmat <- cov(in.df[in.df$OUTL != TRUE, 1:dimension])
    
    print(sprintf("Iter  = %s, cost=%s", i, det(newcovmat)))
    if(!is.null(prevCov) && abs(det(newcovmat) - det(prevCov)) < eps)
    {
      break
    }
    
    prevCov <- covmat
  } #while
    
  m <- as.data.frame(m); colnames(m) <- c("Iterations", "cost_determinant")
  return(list("result" = in.df, "iterations" = m))
}


res      <- findoutliers.mahlanobis(in.df = datf, dimension = dimension, 0.05)
datf.new <- res$result 
ggplot(datf.new, aes(x=x, y=y, color = OUTL)) + geom_point(shape=1)


#--- Let us now try radially shifting the outlier data
hist(datf.new[datf.new$OUTL != TRUE, "mh.dist"])

mh.dist.max <- max(datf.new[datf.new$OUTL != TRUE, "mh.dist"])
center      <- matrix(apply(datf.new[datf.new$OUTL != TRUE, 1:dimension], 2, mean, na.rm = TRUE), ncol = 1)
covmat      <- cov(datf.new[datf.new$OUTL != TRUE, 1:dimension])

bringToBoundary <- function(pointVec, center, covmat, mh.dist.max)
{
  eps <- 1e-10
   mhcurr    <- mahalanobis(x = pointVec, center = center, cov = covmat)

   helper.f <- function(lambda)
   {
     pnew  = lambda* matrix(center, nrow=1) + (1-lambda) * pointVec
     mahalanobis(x = pnew, center = center, cov = covmat) - mh.dist.max
   }

   lambda1  <- 1
   lambda2  <- 0
   lambdaNew <- 1
   
   while(1)
   {
     print(sprintf("lamda1=%s, lambda2=%s, cost=%s",lambda1,lambda2,helper.f(lambdaNew)))
       if(abs(helper.f(lambdaNew)) < eps)
       {
         break
       }
       
       lambdaNew = 0.5 * (lambda1 + lambda2)
       if( helper.f(lambdaNew) * helper.f(lambda1) < 0)
       {
         lambda2  = lambda1
         lambda1  = lambdaNew
       } else if (helper.f(lambdaNew) * helper.f(lambda2) < 0){
         lambda1  = lambda2
         lambda2  = lambdaNew
       } 
     }#while
   
   lambdaNew* matrix(center, nrow=1) + (1-lambdaNew) * pointVec
}
   


head(datf.new)

moved <- data.frame(t(apply(datf.new[datf.new$OUTL == TRUE,1:dimension],1, 
                function(vec){
                              bringToBoundary(pointVec = matrix(vec, nrow = 1), center, covmat, mh.dist.max)
                              })))
names(moved) <- c("x", "y")
moved[["TYPE"]] = "MOVED"
orig.and.moved <- datf.new[datf.new$OUTL == TRUE, 1:dimension]
orig.and.moved[["TYPE"]] = "OLD"
orig.and.moved <- rbind(orig.and.moved,moved )

ggplot(orig.and.moved, aes(x=x, y=y, color = TYPE)) + geom_point(shape=1)
