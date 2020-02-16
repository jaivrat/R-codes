#R Manipulate

library(ggplot2)
library(manipulate)

perimeter = 100

#---------------------------------------------
getResult <- function(len)
{
  x.len <- len
  y.len <- perimeter/2 - x.len
  #tmp <- data.frame(x = 1:len, y = runif(len))
  p   <- ggplot()+
         #ggplot(tmp, aes(x = x, y = y)) + 
         geom_point() + 
         xlim(0, perimeter/2) + ylim(0, perimeter/2) + 
         # Add horizontal line at y = 2O
         geom_hline(yintercept=y.len, linetype="dotted", color = "blue") +
         geom_text(aes(y=y.len, label=sprintf("%s", y.len), x=x.len/2), 
              colour="blue", vjust = -1.2)+
         # Add vertical line at y = 2O
         geom_vline(xintercept =x.len, linetype="dotted", color = "blue")+
         geom_text(aes(x=x.len, label=sprintf("%s", x.len), y=y.len/2), 
                   colour="blue", angle=90, vjust = 1.2 )+
         geom_text(aes(x=perimeter/4, label=sprintf("Area: %s, Peri:%s", x.len * y.len, 2*(x.len + y.len)), y=perimeter/2), 
                   colour="red", angle=0 )+
         geom_rect(mapping=aes(xmin=0, xmax=x.len, ymin=0, ymax=y.len, colour = "blue"), color="black", alpha=0.5) +
         ggtitle("Mr. Shukla(Nutty Prof)")
    
  p
}


#---------------------------------------------
manipulate(
  getResult(len),  
  len = slider(0,perimeter/2, step = 0.01, initial = 10))


#---------------------------------------------
