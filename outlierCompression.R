
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

ggplot(datf, aes(x=x, y=y, color="green")) + geom_point(shape=1) + 
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


res      <- findoutliers.mahlanobis(in.df = datf, dimension = 2, 0.05)
datf.new <- res$result 
ggplot(datf.new, aes(x=x, y=y, color = OUTL)) + geom_point(shape=1)


#









#1. Fist let us generate a data normally distributed with some mean and variance.
#2. We should also generate a t-distributed data

set.seed(1)
N = 10000
good <- rnorm(n = N, mean = 0, sd = 1)
N_bad <- 50
bad  <- rnorm(n = N_bad, mean = 0, sd = 5)
dat <- rbind( data.frame("x" = good, "categ" = "good"),
            data.frame("x" = bad , "categ" = "bad"))


head(dat)

#---------------------- Density Plots --------------------------#
library(ggplot2)
# Density plots with semi-transparent fill
ggplot(dat, aes(x=x, fill=categ)) + geom_density(alpha=.3)

# Density plots with semi-transparent fill
ggplot(dat, aes(x=x)) + geom_density(alpha=.3)


#Create a custom color scale
library(RColorBrewer)
myColors        <- brewer.pal(5,"Set1")
names(myColors) <- levels(dat$categ)
colScale        <- scale_colour_manual(name = "categ",values = myColors)


#---------------------- Scatter Plots 
# Same, but with different colors and add regression lines
ggplot(dat, aes(x=x, y=x, color=categ)) +
  geom_point(shape=1) +
  #scale_colour_hue(l=50) +  # Use a slightly darker palette than normal
  colScale


#---------------------------------------------------------------#
#---------------------------------------------------------------#
#Functions for MAD
DoubleMAD <- function(x, zero.mad.action="warn"){
  # The zero.mad.action determines the action in the event of an MAD of zero.
  # Possible values: "stop", "warn", "na" and "warn and na".
  x         <- x[!is.na(x)]
  m         <- median(x)
  abs.dev   <- abs(x - m)
  left.mad  <- median(abs.dev[x<=m])
  right.mad <- median(abs.dev[x>=m])
  if (left.mad == 0 || right.mad == 0){
    if (zero.mad.action == "stop") stop("MAD is 0")
    if (zero.mad.action %in% c("warn", "warn and na")) warning("MAD is 0")
    if (zero.mad.action %in% c(  "na", "warn and na")){
      if (left.mad  == 0) left.mad  <- NA
      if (right.mad == 0) right.mad <- NA
    }
  }
  return(c(left.mad, right.mad))
}
#- function:
DoubleMADsFromMedian <- function(x, zero.mad.action="warn"){
  # The zero.mad.action determines the action in the event of an MAD of zero.
  # Possible values: "stop", "warn", "na" and "warn and na".
  two.sided.mad <- DoubleMAD(x, zero.mad.action)
  m <- median(x, na.rm=TRUE)
  x.mad <- rep(two.sided.mad[1], length(x))
  x.mad[x > m] <- two.sided.mad[2]
  mad.distance <- abs(x - m) / x.mad
  mad.distance[x==m] <- 0
  return(mad.distance)
}
#check if it works
#pred.tmp[["MAD"]] =  DoubleMADsFromMedian(x) 
#ggplot(data.frame("MAD" = DoubleMADsFromMedian(x)), aes(x=MAD)) + geom_density() + ggtitle(sprintf("%s",datePoints))

#4.
#---- New way of MAD
New_DoubleMADsFromMedian <- function(x, zero.mad.action="warn"){
  # The zero.mad.action determines the action in the event of an MAD of zero.
  # Possible values: "stop", "warn", "na" and "warn and na".
  two.sided.mad <- DoubleMAD(x, zero.mad.action)
  m <- median(x, na.rm=TRUE)
  x.mad <- rep(two.sided.mad[1], length(x))
  x.mad[x > m] <- two.sided.mad[2]
  mad.distance <- (x - m) / x.mad
  mad.distance[x==m] <- 0
  return(mad.distance)
}
#---------------------------------------------------------------#
#---------------------------------------------------------------#




shifterFunc <- function(df, COLUMN)
{
  
  #df will have three cols c("DW_INSTRUMENT_ID", COLUMN , outlflagCol)
  outlflagCol <- "OUTL"
  
  #base case(not outliers)
  if(sum(df[[outlflagCol]] == TRUE) == 0)
  {
    return(df[[COLUMN]])
  }
  
  nonoutl <- df[df[[outlflagCol]] == FALSE, COLUMN]
  s       <- summary(nonoutl)
  q1 <- s[2]; 
  q3 <- s[5];
  med <- s[3];
  
  #....outls...minNonOutl..1Q...med..2Q..MaxNonOutl..outliers
  minNonOutL <- min(df[ df[[outlflagCol]] == FALSE , COLUMN])
  maxNonOutL <- max(df[ df[[outlflagCol]] == FALSE , COLUMN])
  
  minOutL <- min(df[ df[[outlflagCol]] == TRUE , COLUMN], na.rm = TRUE)
  maxOutL <- max(df[ df[[outlflagCol]] == TRUE , COLUMN], na.rm = TRUE)
  
  
  lowdEndIdx  <- (!is.na(df[,COLUMN]) & (df[,COLUMN] < minNonOutL))
  if(sum(lowdEndIdx) > 0)
  {
    lowEndOutls <- df[ lowdEndIdx , COLUMN]
    #(lowEndOutls - minNonOutL)/(minOutL - minNonOutL) = (newValLowEnd - q1)/(maxNonOutL - q1)
    #this means
    newValLowEnd = q1 + ((lowEndOutls - minNonOutL)/(minOutL - minNonOutL)) * (maxNonOutL - q1)
    df[ lowdEndIdx , COLUMN] <- newValLowEnd
  }
  
  
  highEndIdx <- ((!is.na(df[,COLUMN])) & df[,COLUMN] > maxNonOutL)
  if(sum(highEndIdx) > 0)
  {
    highEndOutls <- df[highEndIdx, COLUMN]
    #(highEndOutls - maxNonOutL)/(maxOutL - maxNonOutL) = (newValHighEnd - q3)/(maxNonOutL - q3)
    #this means
    newValHighEnd = q3 + ((highEndOutls - maxNonOutL)/(maxOutL - maxNonOutL)) * (maxNonOutL - q3)
    df[highEndIdx, COLUMN] <- newValHighEnd
  }
  
  return(df[[COLUMN]])
}



absMadOutLier <- function(lb, ub)
{
  function(madDist)
  {
    which(madDist < lb | madDist > ub)
  }
}
dat[["mad"]]  <- New_DoubleMADsFromMedian(dat$x)
flaggingFunction = absMadOutLier(-5,5)
dat[["OUTL"]] <- FALSE
outl.idx <- flaggingFunction(dat[["mad"]])
dat[["OUTL"]][outl.idx] <- TRUE

shifterFunc(dat, COLUMN = "x")





