#https://eight2late.wordpress.com/2017/07/11/a-gentle-introduction-to-logistic-regression-and-lasso-regularisation-using-r/


#load required library
library(mlbench)

#load Pima Indian Diabetes dataset
data("PimaIndiansDiabetes")
head(PimaIndiansDiabetes)
dim(PimaIndiansDiabetes)


#set seed to ensure reproducible results
set.seed(42)

#split into training and test sets
PimaIndiansDiabetes[,"train"] <- ifelse(runif(nrow(PimaIndiansDiabetes))<0.8,1,0)

#separate training and test sets
trainset <- PimaIndiansDiabetes[PimaIndiansDiabetes$train==1,]
testset  <- PimaIndiansDiabetes[PimaIndiansDiabetes$train==0,]

#get column index of train flag
trainColNum <- grep("train",names(trainset))

#remove train flag column from train and test sets
trainset <- trainset[,-trainColNum]
testset  <- testset[,-trainColNum]




#get column index of predicted variable in dataset
typeColNum <- grep("diabetes",names(PimaIndiansDiabetes))

#build model
glm_model <- glm(diabetes~., data = trainset, family = binomial)
summary(glm_model)


#predict probabilities on testset
#type=response gives probabilities, type=class gives class
glm_prob <- predict.glm(glm_model,testset[,-typeColNum],type="response")
#which classes do these probabilities refer to? What are 1 and 0?

#which classes do these probabilities refer to? What are 1 and 0?
contrasts(PimaIndiansDiabetes$diabetes)


#make predictions
##first create vector to hold predictions (we know 0 refers to neg now)
glm_predict <- rep("neg",nrow(testset))
glm_predict[glm_prob>.5] <- "pos"

#confusion matrix
table(pred=glm_predict,true=testset$diabetes)
#glm_predict  neg pos
#neg    90 22
#pos     8 33
#accuracy
mean(glm_predict==testset$diabetes)
#[1] 0.8039216s




#---------------------------------------------------------------------------------
#load required library
library(glmnet)

#convert training data to matrix format
x <- model.matrix(diabetes~.,trainset)

#convert class to numerical variable
y <- ifelse(trainset$diabetes=="pos",1,0)

#perform grid search to find optimal value of lambda
#family= binomial => logistic regression, alpha=1 => lasso
#check docs to explore other type.measure options
cv.out <- cv.glmnet(x, y, alpha=1,family="binomial",type.measure = "mse")
cv.out2<- cv.glmnet(x, y, family="binomial",type.measure = "mse")

#mix of lasso and ridge
#- Log Likelihood function + lambda * [ (1- alpha)*Sum(Beta^2) + alpha * sum(|Beta|)]


#plot result
plot(cv.out)

res.models <- sapply(cv.out$lambda, function(lambda) as.matrix(coef(cv.out,s=lambda)))
row.names <- rownames(coef(cv.out, s=cv.out$lambda.1se))
rownames(res.models) <-row.names
res.models <- rbind(res.models, "loglambda" = log(cv.out$lambda))
res.models <- t(res.models)
res.models <- as.data.frame(res.models[,-2]) #dupliacate Intercept
res.models.melt <- reshape2::melt(res.models, id.vars = "loglambda")
res.models.melt <- res.models.melt[res.models.melt$variable != "(Intercept)", ]
# Map sex to color
ggplot(data=res.models.melt, aes(x=loglambda, y=value, group=variable, colour=variable)) +
  geom_line() +
  geom_vline(xintercept=log(cv.out$lambda.min), aes(color = "red")) +
  geom_vline(xintercept=log(cv.out$lambda.1se), aes(color = "red")) 

coef(cv.out,s=cv.out$lambda.1se)
# 10 x 1 sparse Matrix of class "dgCMatrix"
# 1
# (Intercept) -4.61706681
# (Intercept)  .         
# pregnant     0.03077434
# glucose      0.02314107
# pressure     .         
# triceps      .         
# insulin      .         
# mass         0.02779252
# pedigree     0.20999511
# age     








