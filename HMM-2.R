library(quantmod)

getSymbols("^TWII", src="yahoo", from = "1900-01-01", to = "2015-12-31")

chartSeries(TWII, theme="black")

class(TWII)

x <- as.data.frame(TWII)
head(x)

head(TWII)
#> head(TWII)
# TWII.Open TWII.High TWII.Low TWII.Close TWII.Volume TWII.Adjusted
#1997-07-02   9094.27   9124.30  8988.13    8996.72           0      8996.687
#1997-07-03   9061.60   9061.60  8997.51    9027.63           0      9027.597
#1997-07-04   9144.96   9218.21  9119.25    9192.57           0      9192.536
#> tail(TWII)
#  TWII.Open TWII.High TWII.Low TWII.Close TWII.Volume TWII.Adjusted
# 2015-12-29   8349.27   8351.16  8286.84    8293.91     1186200      8293.878
# 2015-12-30   8313.37   8342.34  8262.52    8279.99     1300600      8279.958
# 2015-12-31   8273.77   8338.06  8258.73    8338.06     1191500      8338.027

TWII_Subset <- window(TWII, start= as.Date("2013-01-01"), end = "2015-12-31")
TWII_Train  <- cbind(TWII_Subset$TWII.Close - TWII_Subset$TWII.Open)
head(TWII_Train)


#Fit HMM though RHmm
library("depmixS4")

data("speed")
class(speed)
head(speed)
set.seed(1)
tst <- runif(4); tst
mod <- depmix(response = rt ~ 1, data = speed, nstates = 2, trstart = tst)
help(depmix)

fm <- fit(mod, emc=em.control(rand=FALSE))
fm
summary(fm)



set.seed(1)
inst = runif(2); inst
mod <- depmix(rt ~ 1, data = speed, 
              nstates = 2, family = gaussian(), 
              transition = ~ scale(Pacc), instart = inst)
fm <- fit(mod, verbose = FALSE, emc=em.control(rand=FALSE))
summary(fm)




#=================================================================
library('depmixS4') #the HMM library we'll use
library('quantmod') 
library(zoo)

setwd('/Users/jvsingh/work/github/rcodes')

#Load out data set (you can download it here), then turn it into a time series object.
#https://inovancetech.com/Data/EURUSD1dHMM.zip

EURUSD1d <- read.csv('EURUSD1d.csv', stringsAsFactors = FALSE, header = TRUE)
Date     <-as.character(EURUSD1d[,1])
DateTS   <- as.POSIXlt(Date, format = "%Y.%m.%d %H:%M:%S") #create date and time objects
TSData   <-data.frame(EURUSD1d[,2:5],row.names=DateTS)
TSData   <-as.xts(TSData) #build our time series data set

ATRindicator <- ATR(TSData[,2:4],n=14) #calculate the indicator
ATR          <- ATRindicator[,2] #grab just the ATR

LogReturns <- log(EURUSD1d$Close) - log(EURUSD1d$Open) #calculate the logarithmic returns

ModelData  <- data.frame(LogReturns,ATR) #create the data frame for our HMM model
ModelData  <- ModelData[-c(1:14),] #remove the data where the indicators are being calculated


#Model
set.seed(1)
         
#We're setting the LogReturns and ATR as our response variables, 
#using the data frame we just built, want to set 3 different regimes, 
#and setting the response distributions to be gaussian.)
HMM <- depmix( list(LogReturns~1, atr~1), 
               data = ModelData,
               nstates=3,
               family=list(gaussian(),gaussian())) 

HMMfit <- fit(HMM, verbose = FALSE) #fit our model to the data set

#we can compare the log Likelihood as well 
#as the AIC and BIC values to help choose our model
print(HMMfit) 

summary(HMMfit)



HMMpost <- depmixS4::posterior(HMMfit) #find the posterior odds for each state over our data set

