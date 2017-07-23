#install.packages("ISLR")

library(ISLR)
data(Smarket)

str(Smarket)

# This data set consists of percentage returns for the S&P 500 stock index over 1, 250 days, 
# from the beginning of 2001 until the end of 2005.

# For each date, we have recorded the percentage returns for each of the five previous 
# trading days, Lag1 through Lag5

# We have also recorded Volume (the number of shares traded on the previous day, in billions), 
# Today (the percentage return on the date in question) and 
# Direction (whether the market was Up or Down on this date).

names(Smarket)
dim(Smarket)
summary(Smarket)

head(Smarket)
pairs(Smarket)

cor(Smarket)
cor(Smarket[,-9])


#----------------------------------------------------------#
#Next, we will fit a logistic regression model in order to predict 
# Direction using Lag1 through Lag5 and Volume

glm.fit = glm(Direction~Lag1+Lag2+Lag3+Lag4+Lag5+Volume , data=Smarket ,family=binomial)
summary(glm.fit)
#The smallest p-value here is associated with Lag1. The negative coefficient 
#for this predictor suggests that if the market had a positive return yesterday, 
#then it is less likely to go up today. However, at a value of 0.15, the p-value 
#is still relatively large, and so there is no clear evidence of a real association 
#between Lag1 and Direction.



glm.probs = predict(glm.fit,type="response")
glm.probs[1:10]

constrasts(Smarket$Direction)

glm.pred = rep("Down", 1250)
glm.pred[glm.probs >.5] = "Up"

tab <- table(glm.pred,Smarket$Direction)
tab

#Accuracy: 0.5216
sum(diag(tab))/sum(tab)


#---------- test/train -----
train=(Smarket$Year <2005)
Smarket.2005= Smarket[! train ,]
dim(Smarket.2005)
Direction.2005=Smarket$Direction[!train]
glm.fit=glm(Direction ~ Lag1+Lag2+Lag3+Lag4+Lag5+Volume , data=Smarket ,family=binomial,subset=train)
glm.probs=predict(glm.fit,Smarket.2005,type="response")
glm.pred=rep("Down",252)
glm.pred[glm.probs >.5]="Up"
table(glm.pred,Direction.2005)
mean(glm.pred==Direction.2005)
mean(glm.pred!=Direction.2005) #The results are rather disappointing: the test error rate is 52%,



#-------------------------------------------------
#We recall that the logistic regression model had very underwhelming p- values associated with all of the predictors, and that the smallest p-value, though not very small, corresponded to Lag1. Perhaps by removing the variables that appear not to be helpful in predicting Direction, we can obtain a more effective model. After all, using predictors that have no relationship with the response tends to cause a deterioration in the test error rate (since such predictors cause an increase in variance without a corresponding decrease in bias), and so removing such predictors may in turn yield an improvement. Below we have refit the logistic regression using just Lag1 and Lag2, which seemed to have the highest predictive power in the original logistic regression model.

glm.fit=glm(Direction ~ Lag1+Lag2,data=Smarket ,family=binomial, subset=train)
glm.probs=predict(glm.fit,Smarket.2005,type="response")
glm.pred=rep("Down",252)
glm.pred[glm.probs >.5]="Up"
table(glm.pred,Direction.2005)

mean(glm.pred==Direction.2005)
#Now the results appear to be more promising: 56 % of the daily movements have been correctly predicted.



#--------------------------------------------------#
#----Linear Discriminant Analysis
library(MASS)
?lda
lda.fit=lda(Direction~Lag1+Lag2,data=Smarket ,subset=train)
lda.fit







