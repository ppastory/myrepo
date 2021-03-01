####################################################################
###########        R code for Simulation methods         ###########
####################################################################

####################################################################
###              Set-up for use                                  ###
####################################################################
rm(list = ls())   # Clear workspace

## Set input file directory

#let's check what is the current directory
getwd()
#it's fiiine !

#We don't have to sent the current directory because it is our directory myrepo
#C:/Users/ppastory/Documents/programming/myrepo
#setwd("C:\\Users\\...") #Windows

## Set output file directory: change to your own directory
#I want the output directory to be the same 
output = getwd() 
#so now everything is in one directory

## Install required packages (code is only executed if the packages are not yet installed)
if(!require(plm)){install.packages("plm")} ## panel data package
#install.packages("expm")


## Load required packages
library(plm) 
#library(expm)

## Load your own functions
source("Case1_Functions.R")

####################################################################
### Monte Carlo simulation for mean of Chi2 distributed variable
####################################################################

#Setting the seeds so that the simulation gives us the same results
set.seed(123)
#n is the size of the sample
T <- 25
#repl number of replication
repl <- 10000
#df is the number of degrees of freedom
df <- 1

#initialising the vector
X_bar <- rep(0,repl)
X_var  <- rep(0,repl)
y_bar <- rep(0,repl)

#initialising the parameters value
b_0 <- 10
b_1 <- 1
beta <- 
sigma2 <- 1
#draw only once the X from a normal distribution
xsim <- rnorm(T,0,1)
#put the constant in X to have right dimension 
X <- as.matrix(cbind(Cnst=1,xsim))
#let's put the beta in matrix form
beta <- as.matrix(rbind(b_0,b_1))
#let's get some errors
e <- rnorm(T,0,sigma2)

#now I can have my y !
Y <- X%*%beta + e

#I am initialising the vectors in which beta and other stuff will arrive
beta_0 <- rep(0,repl)
beta_1 <- rep(0,repl)
stdvs_0 <- rep(0,repl)
stdvs_1 <- rep(0,repl)

beta1_test <- as.matrix(t(c(1,0.95,0.90,0.75,0.5)))

ttest_matrix <- matrix(0,repl, length(beta1_test))

for (j in 1:length(beta_test)) {

  for (i in 1:repl) {
  #stochastic X: X = rnorm(T,0,sigma2)
  #let's get some errors
  e <- rnorm(T,0,sigma2)
  #now I can have my y !
  Y <- X%*%beta + e
  
  OLS_out <- OLS_own(Y,X,0) 
  
  beta_0[i] <- OLS_out$estimation[1,1]
  beta_1[i] <- OLS_out$estimation[2,1]
  
  stdvs_0[i] <- OLS_out$estimation[1,2]/sqrt(T)
  stdvs_1[i] <- OLS_out$estimation[2,2]/sqrt(T)
  
  ttest_matrix[i,j] <- (beta_1[i] - beta1_test[j])/stdvs_1[i]
  print(i,j)
  }

}

colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)


beta_0_bar <- mean(beta_0)
beta_1_bar <- mean(beta_1)
stdvs_0_bar <- mean(stdvs_0)
stdvs_1_bar <- mean(stdvs_1)

#let's do a little table that compares the MC with the population
table <- rbind(cbind(b_0,beta_0_bar),cbind(b_1,beta_1_bar),cbind(0,stdvs_0_bar),cbind(0,stdvs_1_bar))

colnames(table) <- c("population","montecarlo")
table

ttest <- rep(0,length(beta_test))



  
  
  


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
