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
#Case 2.2. MC with known form
####################################################################

#Setting the seeds so that the simulation gives us the same results
set.seed(123)

#we assume we have a ok sample
T <- 200

#repl number of replication
repl <- 10000
#df is the number of degrees of freedom
df <- 1

############################################
######### Initialise Matrix ################

beta_0_mat <- matrix(0,4,5)
colnames(beta_0_mat) <- c("population","estimated beta","std analytical","std numerical","std estimated")


beta_1_mat <- matrix(0,4,5)
colnames(beta_1_mat) <- c("population","estimated beta","std analytical","std numerical","std estimated")


sp_mat <- matrix(0,4,5)
colnames(sp_mat) <- c("size","power B=0.95","power B=0.9","power B=0.75","power B=0.5")



################################################
####Let's start with OLS                  ######
################################################

  
  #initialising the vector
  X_bar <- rep(0,repl)
  X_var  <- rep(0,repl)

  y_bar <- rep(0,repl)
  
  #initialising the parameters value
  b_0 <- 10
  b_1 <- 1
  sigma2 <- 1
  alpa <- 4
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
  beta_0_OLS <- rep(0,repl)
  beta_1_OLS <- rep(0,repl)
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
      
      
      beta_0_0LS[i] <- OLS_out$estimation[1,1]
      beta_1_OLS[i] <- OLS_out$estimation[2,1]
      
      stdvs_0[i] <- OLS_out$estimation[1,2]/sqrt(T)
      stdvs_1[i] <- OLS_out$estimation[2,2]/sqrt(T)
      
      ttest_matrix[i,j] <- (beta_1_OLS[i] - beta1_test[j])/stdvs_1[i] #divided by sqrt T
    }
    
  }
  
  colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)
  
  
  beta_0_bar <- mean(beta_0_0LS)
  beta_1_bar <- mean(beta_1_OLS)
  #let's get the numerical standard errors
  #let's get the numerical standard errors -> truuuue
  
  var_0_num <- var(beta_0_0LS)/T ##divide by T
  var_1_num <- var(beta_1_OLS)/T
  
  stdvs_0_num <- sqrt(var_0_num)
  stdvs_1_num <- sqrt(var_1_num)
  
  #Let get analytical std dev, on a SMALL sample -> truuue
  OLS_std <- OLS_own(Y,X,0)
  stdvs_0_ana <- OLS_std$estimation[1,2]/sqrt(T) 
  stdvs_1_ana <- OLS_std$estimation[2,2]/sqrt(T)
  

  stdvs_0_bar <- mean(stdvs_0)
  stdvs_1_bar <- mean(stdvs_1)
  
 
  #first row of OLS
  table_beta0 <- cbind(b_0,beta_0_bar,stdvs_0_ana,stdvs_0_num,stdvs_0_bar)
  colnames(table_beta0) <- c("population","estimated beta","analytical","numerical","estimated")
  
  #first row of OLS
  table_beta1 <- cbind(b_1,beta_1_bar,stdvs_1_ana,stdvs_1_num,stdvs_1_bar)
  colnames(table_beta1) <- c("population","estimated beta","analytical","numerical","estimated")
  

  
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
  
  
  #Store the result for OLS
  #table with size and power
  sp_mat[1,] <- cbind(size_beta1,t(power_beta1))
  
  #table beta0
  beta_0_mat[1,] <- table_beta0 
  
  
  #table beta1
  beta_1_mat[1,] <-table_beta1

  ################################################
  ####Let's start with OLS White correction ######
  ################################################


  #I am initialising the vectors in which beta and other stuff will arrive
  beta_0_OLSW <- rep(0,repl)
  beta_1_OLSW <- rep(0,repl)
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
      
      OLS_out <- OLS_own(Y,X,1) 
      
      
      beta_0_OLSW[i] <- OLS_out[1,1]
      beta_1_OLSW[i] <- OLS_out[2,1]
      
      stdvs_0[i] <- OLS_out[1,2]/sqrt(T)
      stdvs_1[i] <- OLS_out[2,2]/sqrt(T)
      
      ttest_matrix[i,j] <- (beta_1_OLSW[i] - beta1_test[j])/stdvs_1[i] #divided by sqrt T
    }
    
  }
  
  colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)
  
  
  beta_0_bar <- mean(beta_0_OLSW)
  beta_1_bar <- mean(beta_1_OLSW)

    #let's get the numerical standard errors -> truuuue
  
  var_0_num <- var(beta_0_OLSW)/T ##divide by T
  var_1_num <- var(beta_1_OLSW)/T
  
  stdvs_0_num <- sqrt(var_0_num)
  stdvs_1_num <- sqrt(var_1_num)
  
  #Let get analytical std dev, on a SMALL sample -> truuue
  OLS_std <- OLS_own(Y,X,1)
  stdvs_0_ana <- OLS_std[1,2]/sqrt(T) 
  stdvs_1_ana <- OLS_std[2,2]/sqrt(T)
  
  
  stdvs_0_bar <- mean(stdvs_0)
  stdvs_1_bar <- mean(stdvs_1)
  
  
  #second row of OLSW
  table_beta0 <- cbind(b_0,beta_0_bar,stdvs_0_ana,stdvs_0_num,stdvs_0_bar)
  colnames(table_beta0) <- c("population","estimated beta","analytical","numerical","estimated")
  
  #second row of OLSW
  table_beta1 <- cbind(b_1,beta_1_bar,stdvs_1_ana,stdvs_1_num,stdvs_1_bar)
  colnames(table_beta1) <- c("population","estimated beta","analytical","numerical","estimated")
  
  
  
  #loop over the t-tests and give me the Critical values for each one
  
  #let's do a matrix of critical values
  
  CV_beta1_LB <- quantile(ttest_matrix[,1], c(.025))
  CV_beta1_UB <- quantile(ttest_matrix[,1], c(.975))

  
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
  
  
  #Store the result for OLS with White
  #table with size and power
  sp_mat[2,] <- cbind(size_beta1,t(power_beta1))
  
  #table beta0
  beta_0_mat[2,] <- table_beta0 
  
  
  #table beta1
  beta_1_mat[2,] <-table_beta1
  
  
  
  ################################################
  ####GLS case ######
  ################################################
  
  
  #I am initialising the vectors in which beta and other stuff will arrive
  beta_0_GLS <- rep(0,repl)
  beta_1_GLS <- rep(0,repl)
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
      
      res <- OLS_out$residuals
      
      ### I compute Omega outside of the funciton
      res2       <- res%*%t(res)
      #the non-diagnonal elemets of res2 need to have 0
      #I take away the diagonals
      diagonal <- diag(res2) 
      #I create a matrix of 0 of diam
      res2 <- matrix(0,nrow(res2), ncol(res2)) 
      #I put back the diagonal element in the diagnonals
      #but this time I take the inverse to get P
      diag(res2) <- 1/sqrt(diagonal)
      #res2 is sigma2 omega in the formulas
      P <- res2
      
      omega_1    <- t(P)%*%P
      
      GLS_static = GLS_own (Y,X,omega_1)
      
      
      beta_0_GLS[i] <- GLS_static[1,1]
      beta_1_GLS[i] <- GLS_static[2,1]
      
      stdvs_0[i] <- GLS_static[1,2]/sqrt(T)
      stdvs_1[i] <- GLS_static[2,2]/sqrt(T)
      
      ttest_matrix[i,j] <- (beta_1_GLS[i] - beta1_test[j])/stdvs_1[i] #divided by sqrt T
    }
    
  }
  
  colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)
  
  
  beta_0_bar <- mean(beta_0_GLS)
  beta_1_bar <- mean(beta_1_GLS)
  
  #let's get the numerical standard errors -> truuuue
  
  var_0_num <- var(beta_0_GLS)/T ##divide by T
  var_1_num <- var(beta_1_GLS)/T
  
  stdvs_0_num <- sqrt(var_0_num)
  stdvs_1_num <- sqrt(var_1_num)
  
  #Let get analytical std dev, on a SMALL sample -> truuue
  OLS_std <- OLS_own(Y,X,0)
  res <- OLS_out$residuals
  
  ### I compute Omega outside of the funciton
  res2       <- res%*%t(res)
  #the non-diagnonal elemets of res2 need to have 0
  #I take away the diagonals
  diagonal <- diag(res2) 
  #I create a matrix of 0 of diam
  res2 <- matrix(0,nrow(res2), ncol(res2)) 
  #I put back the diagonal element in the diagnonals
  #but this time I take the inverse to get P
  diag(res2) <- 1/sqrt(diagonal)
  #res2 is sigma2 omega in the formulas
  P <- res2
  
  omega_1    <- t(P)%*%P
  
  GLS_static = GLS_own (Y,X,omega_1)
  
  stdvs_0_ana <- GLS_static[1,2]/sqrt(T)
  stdvs_1_ana <- GLS_static[2,2]/sqrt(T)
  
  #estimated standard errors
  stdvs_0_bar <- mean(stdvs_0)
  stdvs_1_bar <- mean(stdvs_1)
  
  
  #Row of estimation
  table_beta0 <- cbind(b_0,beta_0_bar,stdvs_0_ana,stdvs_0_num,stdvs_0_bar)
  colnames(table_beta0) <- c("population","estimated beta","analytical","numerical","estimated")
  
  #first row of OLS
  table_beta1 <- cbind(b_1,beta_1_bar,stdvs_1_ana,stdvs_1_num,stdvs_1_bar)
  colnames(table_beta1) <- c("population","estimated beta","analytical","numerical","estimated")
  
  
  
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
  
  
  #Store the result for OLS with White
  #table with size and power
  sp_mat[3,] <- cbind(size_beta1,t(power_beta1))
  
  #table beta0
  beta_0_mat[3,] <- table_beta0 
  
  
  #table beta1
  beta_1_mat[3,] <-table_beta1
  
  
  
  ################################################
  ####EGLS case ######
  ################################################
  
  
  #I am initialising the vectors in which beta and other stuff will arrive
  beta_0_EGLS <- rep(0,repl)
  beta_1_EGLS <- rep(0,repl)
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
      
      
      #For EGLS, we assume that Var(miu_{i,t}) is sigma_i^2
      
      n  <- length(Y)
      k  <- ncol(X)

      OLS_out <- OLS_own(Y,X,0) 
      
      res <- OLS_out$residuals
      
      
      x <- X
      y <- log(diag(res%*%t(res)))
      
      ## Run OLS
      xy     <- t(x)%*%y #indeed I need some kind of X_i
      xxi    <- solve(t(x)%*%x)
      coefs  <- as.vector(xxi%*%xy)
      sigma_hat   <- as.vector(x%*%coefs)
      
      sigma_hat <- exp(sigma_hat)
      #now we create the matrix and put sigma_hat as the diagonal of that matrix
      
      omega_hat <- matrix(0,length(res), length(res)) 
      #I put back the diagonal element in the diagnonals
      diag(omega_hat) <- sigma_hat
      
      GLS_static = GLS_own (Y,X,omega_hat)
      
      
      beta_0_EGLS[i] <- GLS_static[1,1]
      beta_1_EGLS[i] <- GLS_static[2,1]
      
      stdvs_0[i] <- GLS_static[1,2]/sqrt(T)
      stdvs_1[i] <- GLS_static[2,2]/sqrt(T)
      
      ttest_matrix[i,j] <- (beta_1_EGLS[i] - beta1_test[j])/stdvs_1[i] #divided by sqrt T
    }
    
  }
  
  colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)
  
  
  beta_0_bar <- mean(beta_0_EGLS)
  beta_1_bar <- mean(beta_1_EGLS)
  
  #let's get the numerical standard errors -> truuuue
  
  var_0_num <- var(beta_0_EGLS)/T ##divide by T
  var_1_num <- var(beta_1_EGLS)/T
  
  stdvs_0_num <- sqrt(var_0_num)
  stdvs_1_num <- sqrt(var_1_num)
  
  #Let get analytical std dev, on a SMALL sample -> truuue
  res <- OLS_own(Y,X,0)$residuals
  
  x <- X
  y <- log(diag(res%*%t(res)))
  
  ## Run OLS
  xy     <- t(x)%*%y #indeed I need some kind of X_i
  xxi    <- solve(t(x)%*%x)
  coefs  <- as.vector(xxi%*%xy)
  sigma_hat   <- as.vector(x%*%coefs)
  
  sigma_hat <- exp(sigma_hat)
  #now we create the matrix and put sigma_hat as the diagonal of that matrix
  
  omega_hat <- matrix(0,length(res), length(res)) 
  #I put back the diagonal element in the diagnonals
  diag(omega_hat) <- sigma_hat
  
  GLS_static = GLS_own (Y,X,omega_hat)
  
  stdvs_0_ana <- GLS_static[1,2]/sqrt(T)
  stdvs_1_ana <- GLS_static[2,2]/sqrt(T)
  
  #estimated standard errors
  stdvs_0_bar <- mean(stdvs_0)
  stdvs_1_bar <- mean(stdvs_1)
  
  
  #Row of estimation
  table_beta0 <- cbind(b_0,beta_0_bar,stdvs_0_ana,stdvs_0_num,stdvs_0_bar)
  colnames(table_beta0) <- c("population","estimated beta","analytical","numerical","estimated")
  
  #first row of OLS
  table_beta1 <- cbind(b_1,beta_1_bar,stdvs_1_ana,stdvs_1_num,stdvs_1_bar)
  colnames(table_beta1) <- c("population","estimated beta","analytical","numerical","estimated")
  
  
  
  #loop over the t-tests and give me the Critical values for each one
  
  #let's do a matrix of critical values
  
  CV_beta1_LB <- quantile(ttest_matrix[,1], c(.025))
  CV_beta1_UB <- quantile(ttest_matrix[,1], c(.975))
  
  
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
  
  
  #Store the result for OLS with White
  #table with size and power
  sp_mat[4,] <- cbind(size_beta1,t(power_beta1))
  
  #table beta0
  beta_0_mat[4,] <- table_beta0 
  
  
  #table beta1
  beta_1_mat[4,] <-table_beta1
  
  
  
  
  
  
  
 