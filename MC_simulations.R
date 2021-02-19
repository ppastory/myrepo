####################################################################
###########        R code for Simulation methods         ###########
####################################################################

####################################################################
###########             Set-up for use                   ###########
####################################################################

rm(list = ls())   # Clear workspace

## For Windows
setwd("C:\\users\\...") # replace the dots by your own path, use getwd() to see current path

## For Mac
## setwd("/users/...") # replace the dots by your own path, use getwd() to see current path

####################################################################
### Monte Carlo simulation for mean of Chi2 distributed variable
####################################################################

set.seed(123)
n = 10
repl = 100000
df = 1

X_bar = rep(0,repl)
X_var  = rep(0,repl)

for (i in 1:repl) {
  X = rchisq(n,df)
  X_bar[i] = mean(X)
  X_var[i] = var(X)/n
}

mean_X_bar = mean(X_bar)
var_X_bar  = var(X_bar)
mean_est_var_X_bar = mean(X_var)
skew_X_bar = mean(((X_bar-mean_X_bar)/sqrt(var_X_bar))^3)
kurt_X_bar = mean(((X_bar-mean_X_bar)/sqrt(var_X_bar))^4)
jb_X_bar = (repl/6)*(skew_X_bar^2 + (1/4)*(kurt_X_bar-3)^2)
p_jb_X_bar = pchisq(jb_X_bar, df=2, lower.tail=FALSE)

####################################################################
### Bootstrap simulation for mean of Chi2 distributed variable
####################################################################

set.seed(123)
n = 10
repl = 100000

data = read.table("data_sim1.csv")
X <- data[1:10,1]

X_bar = rep(0,repl)
X_var  = rep(0,repl)

for (i in 1:repl) {
  sample(X,replace = TRUE)
  X_bar[i] = mean(X)
  X_var[i] = var(X)/n
}

mean_X_bar = mean(X_bar)
var_X_bar  = var(X_bar)
mean_est_var_X_bar = mean(X_var)
skew_X_bar = mean(((X_bar-mean_X_bar)/sqrt(var_X_bar))^3)
kurt_X_bar = mean(((X_bar-mean_X_bar)/sqrt(var_X_bar))^4)
jb_X_bar = (repl/6)*(skew_X_bar^2 + (1/4)*(kurt_X_bar-3)^2)
p_jb_X_bar = pchisq(jb_X_bar, df=2, lower.tail=FALSE)
