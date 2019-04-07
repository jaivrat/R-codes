#install.packages("ldhmm")
library(ldhmm)
library(e1071)

ts = ldhmm.ts_log_rtn ("spx" , on ="days")
#ts$d contains the vector of dates, 
#ts$x contains the vector of log-returns, 
#ts$p contains the vector of prices

DistributionUtils::kurtosis(ts$x)

sapply(0:10, function(drop) e1071::kurtosis(ldhmm.drop_outliers(ts$x , drop)))

#---- Lots of Auto-correlation here. This package needs to address autocorrelation
#The autocorrelation of absolute returns is more than 20%. This feature is 
#called volatility clustering. The markets ACF data can be calculated as
ldhmm.ts_abs_acf(ts$x, drop=0, lag.max =6)
#[1] 0.2444656 0.2737506 0.2502788 0.2434303 0.2809230 0.2389925
ldhmm.ts_abs_acf(ts$x, drop=10, lag.max =6)
#[1] 0.2271080 0.2461340 0.2290556 0.2375244 0.2520413 0.2273132



#Attempt1: Two states by normal distribution: normal and crash regime
#We assume Stationary= TRUE
#As the first experiment, we shall set up the analysis using the normal 
#distribution, which has two parameters: mu and sigma
#(??1, ??1) is approximately (0.0006, 0.01) and (??2, ??2) is approximately (???0.0007, 0.02)
mu_1 = 0.0006; sigma_1 = 0.01
mu_2 = -0.0007; sigma_2 = 0.02
m = 2 #states
gamma0 = ldhmm.gamma_init(m)
param0 = matrix ( c( mu_1 , sigma_1 ,
                     mu_2 , sigma_2 ) , m , 2 , byrow = TRUE )
h <- ldhmm(m, param0, gamma0, stationary=TRUE)
hd <- ldhmm.mle(h , ts$x , decode = TRUE , print.level =2)
#mixing distribution param
hd@param

#Transition probability matrix
hd@gamma

#Stationary state probability vector
hd@delta

hd@states.prob
#However kurtosis is not right
ldhmm.ld_stats(hd)

#Compare this to the statistics of the data classified in each state, as shown below
hd@states.local.stats


#------------------------------------------------------------
#Two states lambda distribution
#------------------------------------------------------------
m = 2
gamma0   <- ldhmm.gamma_init(m)
lambda_1 <- 1; lambda_2 <- 1 #guesses
param0   <- matrix( c( mu_1, sigma_1, lambda_1,
                       mu_2, sigma_2, lambda_2) , m , 3 , byrow = TRUE )
h <- ldhmm(m , param0 , gamma0 , stationary = TRUE )
#Train using MLE
# decode = TRUE : the result object is further 
#                 enriched with the decoded state 
# information based on ts$x
hd <- ldhmm.mle(h, ts$x, decode = TRUE, print.level=2)

#mixing distribution param: 
hd@param
#mu       sigma   lambda
#[1,]  0.0006178873 0.006958494 1.418622
#[2,] -0.0003601983 0.013167710 1.710063
#We see above that both state lambdas > 1 => both state returns are leptokurtic
#Both lambda greater than 1 means leptokurtic

#Transition probability matrix
hd@gamma
#> hd@gamma
#[,1]        [,2]
#[1,] 0.99339008 0.006609922
#[2,] 0.01631603 0.983683971
#> 
  
#Stationary state probability vector
hd@delta
#1] 0.7116838 0.2883162
#According to the stationary state probability vector delta, the market spent 70% of its
#time in normal state.

hd@states.prob

#The statistics of the data classified in each state can be shown by
hd@states.local.stats
#              mean          sd  kurtosis    skewness length
#[1,]  0.0005747770 0.006296764  3.821433 -0.05581264  11977
#[2,] -0.0004506958 0.015358800 17.010688 -0.79468256   4628

#the theoretical statistics of each state is
ldhmm.ld_stats(hd)
#              mean          sd kurtosis
#[1,]  0.0006178873 0.006326603 3.990262
#[2,] -0.0003601983 0.014766337 4.890133
 
#3.990262 ~3.821433 above but second state is very different. This is because very
# large outliers in the data that cannot be captured by two state model.
#We need to drop outliers first, package has the utility to drop outliers
ldhmm.calc_stats_from_obs(hd , drop=11)


ldhmm.oxford_man_plot_obs(hd)
#https://realized.oxford-man.ox.ac.uk/images/oxfordmanrealizedvolatilityindices.zip

#Check by sumulation: we have learnt hmm model, hd above. We sumulate few data from this
hs <- ldhmm.simulate_state_transition ( hd , init =100000)
#measure kurtosis of simulated data
e1071::kurtosis(hs@observations)
#5.624379 : This matches after dropping 10 outliers of actual data=>model is good
#sapply(0:10, function(drop) e1071::kurtosis(ldhmm.drop_outliers(ts$x , drop)))

#The compute-intensive simulation that generates the model's absolute ACF is
#packaged into the ldhmm.simulate_abs_acf function:
ldhmm.simulate_abs_acf(hd , n =100000 , lag.max =1)
#0.1504203
#The ACF number is too low, it is nowhere near 23% that we are expecting. 
#Thus the two-state HMM has some degree of volatility clustering but can not 
#fully reflect the feature for SPX
#We must increase number of states

#We can build more than 3 states models and look for best AIC BIC 
hd@AIC
hd@BIC
head(t(hd@states.prob))
vit <- ldhmm.viterbi(hd, ts$x)
table(vit)
head(vit) # should match with head(t(hd@states.prob)) above



#------------------------------------------------------------
#10 states lambda distribution
#------------------------------------------------------------
m = 10
set.seed(144)
gamma0   <- ldhmm.gamma_init(m)
lambda_s <- sort(runif(m), decreasing = TRUE) + 1  #More kurt to less kurt controlled by lambda
mu_s     <- seq(from = -0.01, to = 0.01, length.out = 10) #neg mean to pos mean
sigma_s  <- seq(from = 0.02, to = 0.001, length.out = 10) #low vol to high vol
#t(sapply(1:m, function(i) c(mu_s[i], sigma_s[i], lambda_s[i])))
#param0   <- matrix( c( mu_1, sigma_1, lambda_1,
#                       mu_2, sigma_2, lambda_2) , m , 3 , byrow = TRUE )
param0 <- t(sapply(1:m, function(i) c(mu_s[i], sigma_s[i], lambda_s[i]))); 
colnames(param0) <- c("mu", "sigma", "lambda")
h.10 <- ldhmm(m , param0 , gamma0 , stationary = TRUE )
#Train using MLE
# decode = TRUE : the result object is further 
#                 enriched with the decoded state 
# information based on ts$x
hd.10 <- ldhmm.mle(h.10, ts$x, decode = TRUE, print.level=2)
#Takes 2.5 hrs to finish
#    user   system  elapsed 
#7691.219 2805.247 9111.364 

#Save the long running model
#saveRDS(object =hd.10, "/Users/jvsingh/work/github/rcodes/hd.10.rds")
#hd.10.saved <- read.saveRDS("/Users/jvsingh/work/github/rcodes/hd.10.rds")
#mixing distribution param: 
hd.10@param

#States may not be necessarily in order:
library(fBasics)
ranks <- rank(apply(hd.10@param, 1, 
           function(e){ fBasics::maxddStats(mean =  e[1], sd = e[2], horizon = 1)}))
#calculate expected maxm drawdown
#7 and #9 look blurred



#mu       sigma   lambda
#[1,]  0.0006178873 0.006958494 1.418622
#[2,] -0.0003601983 0.013167710 1.710063
#We see above that both state lambdas > 1 => both state returns are leptokurtic
#Both lambda greater than 1 means leptokurtic

#Transition probability matrix
floor(hd.10@gamma *10000)/100
#> hd@gamma
#[,1]        [,2]
#[1,] 0.99339008 0.006609922
#[2,] 0.01631603 0.983683971
#> 
  
#Stationary state probability vector
hd@delta
#1] 0.7116838 0.2883162
#According to the stationary state probability vector delta, the market spent 70% of its
#time in normal state.

hd@states.prob

#The statistics of the data classified in each state can be shown by
hd@states.local.stats
#              mean          sd  kurtosis    skewness length
#[1,]  0.0005747770 0.006296764  3.821433 -0.05581264  11977
#[2,] -0.0004506958 0.015358800 17.010688 -0.79468256   4628

#the theoretical statistics of each state is
ldhmm.ld_stats(hd)
#              mean          sd kurtosis
#[1,]  0.0006178873 0.006326603 3.990262
#[2,] -0.0003601983 0.014766337 4.890133
 
#3.990262 ~3.821433 above but second state is very different. This is because very
# large outliers in the data that cannot be captured by two state model.
#We need to drop outliers first, package has the utility to drop outliers
ldhmm.calc_stats_from_obs(hd , drop=11)


ldhmm.oxford_man_plot_obs(hd)
#https://realized.oxford-man.ox.ac.uk/images/oxfordmanrealizedvolatilityindices.zip

#Check by sumulation: we have learnt hmm model, hd above. We sumulate few data from this
hs <- ldhmm.simulate_state_transition ( hd , init =100000)
#measure kurtosis of simulated data
e1071::kurtosis(hs@observations)
#5.624379 : This matches after dropping 10 outliers of actual data=>model is good
#sapply(0:10, function(drop) e1071::kurtosis(ldhmm.drop_outliers(ts$x , drop)))

#The compute-intensive simulation that generates the model's absolute ACF is
#packaged into the ldhmm.simulate_abs_acf function:
ldhmm.simulate_abs_acf(hd , n =100000 , lag.max =1)
#0.1504203
#The ACF number is too low, it is nowhere near 23% that we are expecting. 
#Thus the two-state HMM has some degree of volatility clustering but can not 
#fully reflect the feature for SPX
#We must increase number of states

#We can build more than 3 states models and look for best AIC BIC 
hd@AIC
hd@BIC
head(t(hd@states.prob))
vit <- ldhmm.viterbi(hd, ts$x)
table(vit)
head(vit) # should match with head(t(hd@states.prob)) above

