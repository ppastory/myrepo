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
setwd("C:/Users/ppastory/Documents/programming/myrepo") #Windows

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
#Case 2.1. Monte Carlo simulation for mean of Chi2 distributed variable
####################################################################

#Setting the seeds so that the simulation gives us the same results
set.seed(123)
#n is the size of the sample
sz <- c(25,100,1000,10000,100000)


################################################
####Let's loop over different sample size ######
################################################

#there is one table for beta0

table_final_beta0 <- matrix(0,length(sz),8)
 
colnames(table_final_beta0) <- c("Sample size","population","MC","std_analytical","std_num","std_estimated","Var_ana","Var_estimated")
 
table_final_beta1 <- matrix(0,length(sz), 10)

colnames(table_final_beta1) <- c("Sample size","population","MC","std_analytical","std_num","std_estimated","Var_ana","Var_estimated","2.5%","97.5%")

table_final_sp <- matrix(0,length(sz), 6)

colnames(table_final_sp) <- c("Sample size","size","power B=0.95","power B=0.9","power B=0.75","power B=0.5")

iter <- 0

for (T in sz){
  iter <- iter + 1
  
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

#grid of beta1 
beta1_test <- as.matrix(t(c(1,0.95,0.90,0.75,0.5)))

#ttest matrix with column are different beta1 and 
#the rows are one ttest per simulation

ttest_matrix <- matrix(0,repl, length(beta1_test))

for (j in 1:length(beta1_test)) {

  for (i in 1:repl) {
  #stochastic X: X = rnorm(T,0,sigma2)
  #let's get some errors, we define sigma 2 =1 earlier
  e <- rnorm(T,0,sigma2)
  #now I can have my y !
  Y <- X%*%beta + e
  
  OLS_out <- OLS_own(Y,X,0) 
  
  
  beta_0[i] <- OLS_out$estimation[1,1]
  beta_1[i] <- OLS_out$estimation[2,1]
  
  stdvs_0[i] <- OLS_out$estimation[1,2]/sqrt(T)
  stdvs_1[i] <- OLS_out$estimation[2,2]/sqrt(T)
  
  ttest_matrix[i,j] <- (beta_1[i] - beta1_test[j])/stdvs_1[i] #divided by sqrt T
  }

}

colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)


beta_0_bar <- mean(beta_0)
beta_1_bar <- mean(beta_1)

#let's get the numerical standard errors
stdvs_0_num <- sqrt(var(beta_0))/sqrt(T)
stdvs_1_num <- sqrt(var(beta_1))/sqrt(T)

var_0_num <- var(beta_0)/T ##divide by T?
var_1_num <- var(beta_1)/T


#Let get analytical std dev, on a SMALL sample
OLS_std <- OLS_own(Y,X,0)
stdvs_0_ana <- OLS_std$estimation[1,2]/sqrt(T) 
stdvs_1_ana <- OLS_std$estimation[2,2]/sqrt(T)


x<-X
xxi    <- solve(t(x)%*%x) #this is (X' X)^(-1)
var_01_ana  <- sigma2*diag(xxi)/T
var_0_ana <- var_01_ana[1]
var_1_ana <- var_01_ana[2]

stdvs_0_bar <- mean(stdvs_0)
stdvs_1_bar <- mean(stdvs_1)

var_0_bar <- mean(stdvs_0)^2
var_1_bar <- mean(stdvs_1)^2



#let's see what are the true standard standars
table_std <- rbind(cbind(stdvs_0_ana,stdvs_0_num,stdvs_0_bar),cbind(stdvs_1_ana,stdvs_1_num,stdvs_1_bar))
colnames(table_std) <- c("analytical","numerical","estimated")
table_std

#let's see what are the true variance
table_var <- rbind(cbind(var_0_ana,var_0_num,var_0_bar),cbind(var_1_ana,var_1_num,var_1_bar))
colnames(table_var) <- c("analytical","numerical","estimated")
table_var


#let's do a little table that compares the MC with the population
table_beta <- cbind(b_0,beta_0_bar,b_1,beta_1_bar)

colnames(table_beta) <- c("Beta_0_pop","Beta_0_MC","Beta_1_pop","Beta_1_MC")
table_beta <- cbind(T,table_beta)


#loop over the t-tests and give me the Critical values for each one

#let's do a matrix of critical values

CV_beta1_LB <- quantile(ttest_matrix[,1], c(.025))
CV_beta1_UB <- quantile(ttest_matrix[,1], c(.975))

#

rej_function <- function(ttest,LB,UB){
  if (ttest <= LB || ttest >= UB){
    return(1)
  }else{
    return(0)
  }
}

#initialise matrix of rejection 
rej_matrix <- matrix(0,nrow(ttest_matrix),ncol(ttest_matrix))
#Lets store the 1 and 0 of rejection in a matrix
for (j in 1:5){
  rej_matrix[,j]<-sapply(ttest_matrix[,j],rej_function, LB= CV_beta1_LB,UB= CV_beta1_UB)
}
#this is the size, the mean of column 1 for which beta =1
size_beta1 <- mean(rej_matrix[,1])
#size_beta1

#the power is P(non reject if Beta != 1) -> 1 - P(reject)
power_beta1 <- 1 -colMeans(rej_matrix[,2:5])
#power_beta1

#this is the table in which you have all beta_0
table_final_beta0[iter,] <- cbind(T,b_0,beta_0_bar,stdvs_0_ana,stdvs_0_num,stdvs_0_bar,var_0_ana,var_0_bar)

#this is the table in which you have all beta_1
table_final_beta1[iter,] <- cbind(T,b_1,beta_1_bar,stdvs_1_ana,stdvs_1_num,stdvs_1_bar,var_1_ana,var_1_bar,CV_beta1_LB,CV_beta1_UB)

#table with size and power
table_final_sp[iter,] <- cbind(T,size_beta1,t(power_beta1))

}








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
