###################################################################
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
#Case 2.3. MC with unknown form
####################################################################

#Setting the seeds so that the simulation gives us the same results
set.seed(123)

#we assume we have a ok sample
T <- 500

#repl number of replication
repl <- 1500 #less number of replication to work on the code


############################################
######### Initialise Matrix ################


sp_mat <- matrix(0,2,5)
colnames(sp_mat) <- c("size","power B=0.95","power B=0.9","power B=0.75","power B=0.5")



################################################
####OLS Bootstrap only using one MC draw  ######
################################################


#initialising the vector
X_bar <- rep(0,repl)
X_var  <- rep(0,repl)

y_bar <- rep(0,repl)

#initialising the parameters value
b_0 <- 10
b_1 <- 1
sigma2 <- 1
alpha <- 4
#draw only once the X from a normal distribution
xsim <- rnorm(T,0,1)
#put the constant in X to have right dimension 
X <- as.matrix(cbind(Cnst=1,xsim))
#let's put the beta in matrix form
beta <- as.matrix(rbind(b_0,b_1))

#let's get some errors: they are heteroskedastic
diagonal <- xsim^alpha
#the errors are heteroskedastic
e <- rnorm(T,0,sd = sqrt(diagonal))


#now I can have my y !
Y <- X%*%beta + e

#we need the estimated residuals
OLS_out <- OLS_own(Y,X,0)
res <- OLS_out$residuals

  
#bootstrap replication
brepl <- 100

#I am initialising the vectors in which beta and other stuff will arrive -Pairwise
beta_0_0LSbp <- rep(0,brepl)
beta_1_OLSbp <- rep(0,brepl)

#I am initialising the vectors in which beta and other stuff will arrive- Wild
beta_0_0LSbw <- rep(0,brepl)
beta_1_OLSbw <- rep(0,brepl)

#grid of beta1 
beta1_test <- as.matrix(t(c(1,0.95,0.90,0.75,0.5)))

#ttest matrix with column are different beta1 and 
#the rows are one ttest per simulation

ttest_matrixbp <- matrix(0,repl, length(beta1_test))
ttest_matrixbw <- matrix(0,repl, length(beta1_test))


data <- cbind(Y,X)

for (j in 1:length(beta1_test)) {
  
  
  for (b in 1:brepl) {

    #1 pair bootstrap 
    #I am taking our a boostrap pair of Y and X
    boot_pair = data[unlist(sample(as.data.frame(matrix(1:nrow(data),nrow = 2)),100,replace=T)),]
    #by taking a pair of Y and x -> then x is not deterministic anymore
    
    #My y in the pair is the last column
    Ybp<-as.matrix(boot_pair[,1])
    #My x in the pair is everything but the last column
    Xbp<-boot_pair[,-1]
    
    #let's compute OLS with our pair-wise sample
    OLS_out <- OLS_own(Ybp,Xbp,0) 
    
    beta_0_0LSbp[b] <- OLS_out$estimation[1,1]
    beta_1_OLSbp[b] <- OLS_out$estimation[2,1]
    
    #2 Wild bootstrap 
    #I generate the sample of shocks to get some stochasticity
    #I generate 1 and -1
    shock <- sample(c(-1,1), replace=TRUE, size=T)
    
    #I put a shock to the residuals of the original data
    e_b <- sample(res, replace=TRUE, size=T)
    #let's get the errors: I need element-wise multiplication
    errors_b <- shock*e_b
    
    X_b <- X
      
    #as.matrix(cbind(Cnst=1,sample(X[,2], replace=TRUE, size=200)))
    
    #now I can have my y !
    Y_b <- X%*%beta + errors_b
    
    OLS_out <- OLS_own(Y_b,X_b,0) 
    
    beta_0_0LSbw[b] <- OLS_out$estimation[1,1]
    beta_1_OLSbw[b] <- OLS_out$estimation[2,1]
    
    
  }
  
}

colnames(ttest_matrixbp) <- c(1,0.95,0.90,0.75,0.5)
colnames(ttest_matrixbw) <- c(1,0.95,0.90,0.75,0.5)


#Do the pairwise standard errors

beta_0_barp <- mean(beta_0_0LSbp)
beta_1_barp <- mean(beta_1_OLSbp)

#let's get the numerical standard errors -> true ones
var_0_nump <- var(beta_0_0LSbp) 
var_1_nump <- var(beta_1_OLSbp)

stdvs_0_nump <- sqrt(var_0_nump)
stdvs_1_nump <- sqrt(var_1_nump)



###Do the wild standard errors 

beta_0_barw <- mean(beta_0_0LSbw)
beta_1_barw <- mean(beta_1_OLSbw)

#let's get the numerical standard errors 

var_0_numw <- var(beta_0_0LSbw)
var_1_numw <- var(beta_1_OLSbw)

stdvs_0_numw <- sqrt(var_0_numw)
stdvs_1_numw <- sqrt(var_1_numw)

#first row of OLS
table_std_num <- cbind(stdvs_0_nump,stdvs_1_nump,stdvs_0_numw,stdvs_1_numw)
colnames(table_std_num) <- c("pair std beta0","pair std beta1","wild std beta0","wild std beta1")
table_std_num



################################################
####OLS Bootstrap only using many MC draw  ######
################################################

beta_0_OLS <- rep(0,repl)
beta_1_OLS <- rep(0,repl)

stdvs_0_OLS <- rep(0,repl)
stdvs_1_OLS <- rep(0,repl)

#For every MC simulation, I get an estimated beta and its standard errors !

beta_0_barp<- rep(0,repl) 
beta_1_barp<- rep(0,repl)

stdvs_0bp<- rep(0,repl)
stdvs_1bp<- rep(0,repl)

beta_0_barw<- rep(0,repl)
beta_1_barw<- rep(0,repl)

stdvs_0bw<- rep(0,repl)
stdvs_1bw<- rep(0,repl)

for (j in 1:length(beta1_test)) {
  #for each MC simulation -> I draw an X, a Y an error
  for (i in 1:repl) {

    #let's get some errors: they are heteroskedastic
    diagonal <- xsim^alpha
    #the errors are heteroskedastic
    e <- rnorm(T,0,sd = sqrt(diagonal))    
    #now I can have my y !
    Y <- X%*%beta + e
    
    #we need the estimated residuals
    OLS_out <- OLS_own(Y,X,0)
    res <- OLS_out$residuals
    
    
    print(i)
    
      for (b in 1:brepl) {
        
        #1 pair bootstrap 
        #I am taking our a boostrap pair of Y and X
        boot_pair = data[unlist(sample(as.data.frame(matrix(1:nrow(data),nrow = 2)),100,replace=T)),]
        
        #My y in the pair is the last column
        Ybp<-as.matrix(boot_pair[,1])
        #My x in the pair is everything but the last column
        Xbp<-boot_pair[,-1]
        
        #let's compute OLS with our pair-wise sample
        OLS_out <- OLS_own(Ybp,Xbp,0) 
        
        beta_0_0LSbp[b] <- OLS_out$estimation[1,1]
        beta_1_OLSbp[b] <- OLS_out$estimation[2,1]
        
        
        #2 Wild bootstrap 
        #I generate the sample of shocks to get some stochasticity
        #I generate 1 and -1
        shock <- sample(c(-1,1), replace=TRUE, size=T)
        
        #I put a shock to the residuals of the original data
        e_b <- sample(res, replace=TRUE, size=T)
        #let's get the errors: I need element-wise multiplication
        errors_b <- shock*e_b
        
        X_b <- X
        
        #as.matrix(cbind(Cnst=1,sample(X[,2], replace=TRUE, size=200)))
        
        #now I can have my y !
        Y_b <- X%*%beta + errors_b
        
        OLS_out <- OLS_own(Y_b,X_b,0) 
        
        beta_0_0LSbw[b] <- OLS_out$estimation[1,1]
        beta_1_OLSbw[b] <- OLS_out$estimation[2,1]
        
        
      }
  
    #Do the pairwise standard errors
    #let's get the numerical standard errors -> true ones
    var_0_nump <- var(beta_0_0LSbp) 
    var_1_nump <- var(beta_1_OLSbp)
    
    stdvs_0_nump <- sqrt(var_0_nump)
    stdvs_1_nump <- sqrt(var_1_nump)
    
    #I get a bootstrap beta  
    beta_0_barp[i] <- mean(beta_0_0LSbp)
    beta_1_barp[i] <- mean(beta_1_OLSbp)
    
    #the mean of this standard-errors from pair-sample will be the estimated standard errors
    #This needs to be done manually because of sl 37
    #the function takes into account the degree correction which inflates the variance
    #-> it is not necessary to do it manually
    stdvs_0bp[i] <- stdvs_0_nump
    stdvs_1bp[i] <- stdvs_1_nump
    #Compute the t test with these standard errors
    ttest_matrixbp[i,j] <- (beta_1_barp[i] - beta1_test[j])/stdvs_1bp[i] 
    
    
    ###Do the wild standard errors
    #let's get the numerical standard errors 
    
    var_0_numw <- var(beta_0_0LSbw)
    var_1_numw <- var(beta_1_OLSbw)
    
    stdvs_0_numw <- sqrt(var_0_numw)
    stdvs_1_numw <- sqrt(var_1_numw)
    
    #For every Monte-carlo I get a Wild Bootstrap std errors
    stdvs_0bw[i] <- stdvs_0_numw
    stdvs_1bw[i] <- stdvs_1_numw
    
    #I get a bootstrap beta  
    beta_0_barw[i] <- mean(beta_0_0LSbw)
    beta_1_barw[i] <- mean(beta_1_OLSbw) 
    
    
    #And I get a t stat based on the bootstrap
    ttest_matrixbw[i,j] <- (beta_1_barw[i] - beta1_test[j])/stdvs_1bw[i] 
    
    
    #I need the MC estimated std errors
    OLS_out <- OLS_own(Y,X,0)
    
    beta_0_OLS[i] <- OLS_out$estimation[1,1]
    beta_1_OLS[i] <- OLS_out$estimation[2,1]
    
    #the mean of this standard-errors from pair-sample will be the estimated standard errors
    #This needs to be done manually because of sl 37
    #the function takes into account the degree correction which inflates the variance
    #-> it is not necessary to do it manually
    stdvs_0_OLS[i] <- OLS_out$estimation[1,2]
    stdvs_1_OLS[i] <- OLS_out$estimation[2,2]
    
    
    }   
}

colnames(ttest_matrixbp) <- c(1,0.95,0.90,0.75,0.5)
colnames(ttest_matrixbw) <- c(1,0.95,0.90,0.75,0.5)

#We want to check that the estimated standard error is unbiased
#A. Estimated standard error#
####

stdvs_0_est_MC <- mean(stdvs_0_OLS)
stdvs_1_est_MC <- mean(stdvs_1_OLS)

stdvs_0_est_BP <- mean(stdvs_0bp)
stdvs_1_est_BP <- mean(stdvs_1bp)

stdvs_0_est_BW <- mean(stdvs_0bw)
stdvs_1_est_BW <- mean(stdvs_1bw)

####
#B. True standard errors
####

###
#B.1 Use an analytical formula
###

##Analytical OLS with heteroskedasticity

sigma2 <- 1
alpha <- 4

diagonal <- xsim^alpha
#sigma_omega is the asymptotic variance of the OLS estimator
sigma_omega <- matrix(0,T,T) 
#I put back the diagonal element in the diagnonal  
diag(sigma_omega) <- diagonal
x <- X
xxi    <- solve(t(x)%*%x) #this is (X' X)^(-1)
cov_OLS <- xxi %*% t(x) %*% sigma_omega %*% x %*% xxi

stdvs_0_ana_OLS <- sqrt(cov_OLS[1,1])
stdvs_1_ana_OLS <- sqrt(cov_OLS[2,2])


##Analytical OLSW
OLS_std <- OLS_own(Y,X,1)
stdvs_0_ana_OLSW <- OLS_std[1,2]
stdvs_1_ana_OLSW <- OLS_std[2,2]

###
#B.2 Numerical Beta (MS standard errors)
###

var_0_num <- var(beta_0_OLS)
var_1_num <- var(beta_1_OLS)

stdvs_0_numMC_OLS <- sqrt(var_0_num)
stdvs_1_numMC_OLS <- sqrt(var_1_num)




###store the betas

#beta_0 table
table_beta0 <- cbind(b_0,stdvs_0_est_MC,stdvs_0_numMC_OLS,stdvs_0_est_BP,stdvs_0_est_BW,stdvs_0_ana_OLS,stdvs_0_ana_OLSW)
colnames(table_beta0) <- c("population","estimated std MC","numerical MC std","estimated bootstrap pairwise","estimated wild bootstrap","analytical OLS","analytical OLSW")

#beta_1 table
table_beta1 <- cbind(b_1,stdvs_1_est_MC,stdvs_1_numMC_OLS,stdvs_1_est_BP,stdvs_1_est_BW,stdvs_1_ana_OLS,stdvs_1_ana_OLSW)
colnames(table_beta0) <- c("population","estimated std MC","numerical MC std","estimated bootstrap pairwise","estimated wild bootstrap","analytical OLS","analytical OLSW")



#Bootstap pairwise, size and power !!

#let's do a matrix of critical values

CV_beta1_LBbp <- quantile(ttest_matrixbp[,1], c(.025))
CV_beta1_UBbp <- quantile(ttest_matrixbp[,1], c(.975))

#

rej_function <- function(ttest,LB,UB){
  if (ttest <= LB || ttest >= UB){
    return(1)
  }else{
    return(0)
  }
}

#initialise matrix of rejection 
rej_matrixbp <- matrix(0,nrow(ttest_matrixbp),ncol(ttest_matrixbp))
#Lets store the 1 and 0 of rejection in a matrix
for (j in 1:5){
  rej_matrixbp[,j]<-sapply(ttest_matrixbp[,j],rej_function, LB= CV_beta1_LBbp,UB= CV_beta1_UBbp)
}
#this is the size, the mean of column 1 for which beta =1
size_beta1bp <- mean(rej_matrixbp[,1])
#size_beta1

#the power is P(non reject if Beta != 1) -> 1 - P(reject)
power_beta1bp <- colMeans(rej_matrixbp[,2:5])
#power_beta1


#Store the result for OLS
#table with size and power
sp_mat[1,] <- cbind(size_beta1bp,t(power_beta1bp))


#Bootstrap wild size and power

#loop over the t-tests and give me the Critical values for each one

#let's do a matrix of critical values

CV_beta1_LBbw <- quantile(ttest_matrixbw[,1], c(.025))
CV_beta1_UBbw <- quantile(ttest_matrixbw[,1], c(.975))

#

rej_function <- function(ttest,LB,UB){
  if (ttest <= LB || ttest >= UB){
    return(1)
  }else{
    return(0)
  }
}

#initialise matrix of rejection 
rej_matrixbw <- matrix(0,nrow(ttest_matrixbw),ncol(ttest_matrixbw))
#Lets store the 1 and 0 of rejection in a matrix
for (j in 1:5){
  rej_matrixbw[,j]<-sapply(ttest_matrixbw[,j],rej_function, LB= CV_beta1_LBbw,UB= CV_beta1_UBbw)
}
#this is the size, the mean of column 1 for which beta =1
size_beta1bw <- mean(rej_matrixbw[,1])
#size_beta1

#the power is P(non reject if Beta != 1) -> 1 - P(reject)
power_beta1bw <- colMeans(rej_matrixbw[,2:5])
#power_beta1


#Store the result for OLS
#table with size and power
sp_mat[2,] <- cbind(size_beta1bw,t(power_beta1bw))

