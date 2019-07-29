rm(list = ls())
setwd('/Users/jvsingh/work/github/R-codes')

datf <- read.csv(file = "data/custer_points.csv", header = TRUE, stringsAsFactors = FALSE)

library(ggplot2)



#SVM demo
#Random init of points
x.max <- max(datf$X)
x.min <- min(datf$X)
y.max <- max(datf$Y)
y.min <- min(datf$Y)

#random points
#set.seed(144)
p1 <- c(x.min + runif(1)* (x.max - x.min), y.min + runif(1)*( y.max - y.min))
p2 <- c(x.min + runif(1)* (x.max - x.min), y.min + runif(1)*( y.max - y.min))

sleep_time = 5

for(i in 1:10)
{
  print(i); print("...")
  if(i==1)
  {
   g <-  ggplot(datf, aes(x = X, y = Y, color = category) )+
      geom_point()+
      annotate("point", x = p1[1], y = p1[2], colour = "red",  shape = 3 )+
      annotate("point", x = p2[1], y = p2[2], colour = "blue",  shape = 3 )
   
   print(g)
    Sys.sleep(sleep_time)  
  }
  
  #assignment step
  dists <- t(apply(datf[,1:2], 1, 
                 function(this.row){
                   x = this.row[1]; y <- this.row[2]
                   c(sqrt(sum((this.row[1:2] - p1)^2)), sqrt(sum((this.row[1:2] - p2)^2)))
                 }))
  colours <- ifelse(apply(dists, 1, function(dis) which(min(dis) == dis)) == 1, "red", "blue")
  datf$category <- colours
  g <- ggplot(datf, aes(x = X, y = Y, color = category) )+
    geom_point(aes(colour = category))+scale_colour_identity()+
    annotate("point", x = p1[1], y = p1[2], colour = "red",  shape = 3 )+
    annotate("point", x = p2[1], y = p2[2], colour = "blue",  shape = 3 )
  print(g)
  Sys.sleep(sleep_time)  
  
  
  #Centroid step
  p1 <- c( mean(datf[datf$category == "red", "X"]), mean(datf[datf$category == "red", "Y"]))
  p2 <- c( mean(datf[datf$category == "blue", "X"]), mean(datf[datf$category == "blue", "Y"]))

  g <- ggplot(datf, aes(x = X, y = Y, color = category) )+
    geom_point(aes(colour = category))+scale_colour_identity()+
    annotate("point", x = p1[1], y = p1[2], colour = "red",  shape = 3 )+
    annotate("point", x = p2[1], y = p2[2], colour = "blue",  shape = 3 )
  print(g)
  Sys.sleep(0)
}




