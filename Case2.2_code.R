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

#we assume we have sample of reasonable size
T <- 500

#repl number of replication
repl <- 1000 #less number of replication to work on the code


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
  
  #draw only once the X from a normal distribution
  #X is fixed over repeated samples
  xsim <- rnorm(T,0,1)
  #put the constant in X to have right dimension 
  X <- as.matrix(cbind(Cnst=1,xsim))
  #let's put the beta in matrix form
  beta <- as.matrix(rbind(b_0,b_1))
  
  #let's get some errors
  #sigma2 is 1
  sigma2 <- 1
  alpha <- 4
  #create the vector of errors
  e <- rnorm(T,0,sigma2*xsim^alpha)
  
  #now I can have my y !
  Y <- X%*%beta + e

  #I am initialising the vectors in which beta and other stuff will arrive
  beta_0_0LS <- rep(0,repl)
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
      e <- rnorm(T,0,sigma2*xsim^alpha)
      #now I can have my y !
      Y <- X%*%beta + e
      
      OLS_out <- OLS_own(Y,X,0) 
      
      
      beta_0_0LS[i] <- OLS_out$estimation[1,1]
      beta_1_OLS[i] <- OLS_out$estimation[2,1]
    
      
      stdvs_0[i] <- OLS_out$estimation[1,2] #this is se(B^MC) sl 35
      stdvs_1[i] <- OLS_out$estimation[2,2]
      
      ttest_matrix[i,j] <- (beta_1_OLS[i] - beta1_test[j])/stdvs_1[i] 
    }
    print(mean(stdvs_0))
  }
  
  colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)
  
  
  beta_0_est <- mean(beta_0_0LS)
  beta_1_est <- mean(beta_1_OLS)

  #Numerical Standard errors
  
  var_0_num <- var(beta_0_0LS)
  var_1_num <- var(beta_1_OLS)
  
  stdvs_0_num <- sqrt(var_0_num)
  stdvs_1_num <- sqrt(var_1_num)
  
  #Analytical standard errors
  
  sigma2_i <- sigma2*xsim^alpha
  #sigma2_i is the diagonal of sigma2 omega
  sigma_omega <- matrix(0,T,T) 
  #I put back the diagonal element in the diagnonals
  diag(sigma_omega) <- sigma2_i
  #sigma_omega is the asymptotic variance of the OLS estimator
  
  #using expression in slide 53
  x <- X
  xxi    <- solve(t(x)%*%x) #this is (X' X)^(-1)
  cov_OLS <- xxi %*% t(x) %*% sigma_omega %*% x %*% xxi
  
  stdvs_0_ana <- sqrt(cov_OLS[1,1])
  stdvs_1_ana <- sqrt(cov_OLS[2,2])
  

  stdvs_0_est <- mean(stdvs_0)
  stdvs_1_est <- mean(stdvs_1)
  
 
  #first row of OLS
  table_beta0 <- cbind(b_0,beta_0_est,stdvs_0_ana,stdvs_0_num,stdvs_0_est)
  colnames(table_beta0) <- c("population","estimated beta","analytical","numerical","estimated")
  
  #first row of OLS
  table_beta1 <- cbind(b_1,beta_1_est,stdvs_1_ana,stdvs_1_num,stdvs_1_est)
  colnames(table_beta1) <- c("population","estimated beta","analytical","numerical","estimated")
  

  
  ###
  #C Compute size and power
  ###
  
  #loop over the t-tests and give me the Critical values for each one
  colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)
  
  #let's do a matrix of critical values, alpha is 5%
  #the t statistics follows a student t with 2 degrees of freedom 
  
  CV_beta1 <- qt(p=.05/2, df=T-2, lower.tail=FALSE)
  
  rej_function <- function(ttest,cv){
    if (abs(ttest) > cv){
      return(1)
    }else{
      return(0)
    }
  }
  
  
  #initialise matrix of rejection 
  rej_matrix <- matrix(0,nrow(ttest_matrix),ncol(ttest_matrix))
  #Lets store the 1 and 0 of rejection in a matrix
  for (j in 1:5){
    rej_matrix[,j]<-sapply(ttest_matrix[,j],rej_function, cv= CV_beta1)
  }
  #this is the size, the mean of column 1 for which beta = 1
  size_beta1 <- mean(rej_matrix[,1])
  #there is a size bias if the critical values are wrong !
  
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
  
  #here OLS standard errors will be biased and inconsistent

  OLS_out <- OLS_own(Y,X,1) 
  
  
  beta_0_OLSW[i] <- OLS_out[1,1]
  
  
  
  ################################################
  ####Let's start with OLS White correction ######
  ################################################
  sigma2 <- 1

  #I am initialising the vectors in which beta and other stuff will arrive
  beta_0_OLSW <- rep(0,repl)
  beta_1_OLSW <- rep(0,repl)
  stdvs_0 <- rep(0,repl)
  stdvs_1 <- rep(0,repl)
  sigma2 <- 1
  
  #grid of beta1 
  beta1_test <- as.matrix(t(c(1,0.95,0.90,0.75,0.5)))
  
  #ttest matrix with column are different beta1 and 
  #the rows are one ttest per simulation
  
  ttest_matrix <- matrix(0,repl, length(beta1_test))
  
 # for (j in 1:length(beta1_test)) {
    
    for (i in 1:repl) {
      #stochastic X: X = rnorm(T,0,sigma2)
      #let's get some errors, we define sigma 2 =1 earlier
      e <- rnorm(T,0,sigma2*xsim^alpha)
      #now I can have my y !
      Y <- X%*%beta + e
      
      OLS_out <- OLS_own(Y,X,1) 
      
      
      beta_0_OLSW[i] <- OLS_out[1,1]
      beta_1_OLSW[i] <- OLS_out[2,1]
      
      
      stdvs_0[i] <- OLS_out[1,2]
      stdvs_1[i] <- OLS_out[2,2]
      
      ttest_matrix[i,j] <- (beta_1_OLSW[i] - beta1_test[j])/stdvs_1[i] 
    }
    
 # }
  
  colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)
  
  
  beta_0_est <- mean(beta_0_OLSW)
  beta_1_est <- mean(beta_1_OLSW)

  #Numerical Standard errors
  
  var_0_num <- var(beta_0_OLSW) 
  var_1_num <- var(beta_1_OLSW)
  
  stdvs_0_num <- sqrt(var_0_num)
  stdvs_1_num <- sqrt(var_1_num)
  
  #Analytical Standard errors #careful here I don't assum sigma2 =1 !!
  
  OLS_std <- OLS_own(Y,X,1)
  stdvs_0_ana <- OLS_std[1,2]
  stdvs_1_ana <- OLS_std[2,2]
  
  
  stdvs_0_est <- mean(stdvs_0)
  stdvs_1_est <- mean(stdvs_1)
  
  
  #second row of OLSW
  table_beta0 <- cbind(b_0,beta_0_est,stdvs_0_ana,stdvs_0_num,stdvs_0_est)
  colnames(table_beta0) <- c("population","estimated beta","analytical","numerical","estimated")
  
  #second row of OLSW
  table_beta1 <- cbind(b_1,beta_1_est,stdvs_1_ana,stdvs_1_num,stdvs_1_est)
  colnames(table_beta1) <- c("population","estimated beta","analytical","numerical","estimated")
  
  
  
  ###
  #C Compute size and power
  ###
  
  #loop over the t-tests and give me the Critical values for each one
  colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)
  
  #let's do a matrix of critical values, alpha is 5%
  #the t statistics follows a student t with 2 degrees of freedom 
  
  CV_beta1 <- qt(p=.05/2, df=T-2, lower.tail=FALSE)
  
  rej_function <- function(ttest,cv){
    if (abs(ttest) > cv){
      return(1)
    }else{
      return(0)
    }
  }
  
  
  #initialise matrix of rejection 
  rej_matrix <- matrix(0,nrow(ttest_matrix),ncol(ttest_matrix))
  #Lets store the 1 and 0 of rejection in a matrix
  for (j in 1:5){
    rej_matrix[,j]<-sapply(ttest_matrix[,j],rej_function, cv= CV_beta1)
  }
  #this is the size, the mean of column 1 for which beta = 1
  size_beta1 <- mean(rej_matrix[,1])
  #there is a size bias if the critical values are wrong !
  
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
  
  #Here the white correction will correct for heteroskedasticty 
  #standard errors will be biased (in small sample) but consistent
  #it is inefficient -> the GLS standard errors will be smaller
  
  
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
  
 # for (j in 1:length(beta1_test)) {
    
    for (i in 1:repl) {
      sigma2 <-1
      #stochastic X: X = rnorm(T,0,sigma2)
      #let's get some errors, we define sigma 2 =1 earlier
      e <- rnorm(T,0,sigma2*xsim^alpha)
      #now I can have my y !
      Y <- X%*%beta + e
      
      OLS_out <- OLS_own(Y,X,0) 
      
      res <- OLS_out$residuals
      
      ### I observe the variance!
      res2       <- res%*%t(res)
      #the non-diagnonal elemets of res2 need to have 0
      #I take away the diagonals
      #diagonal <- diag(res2) 
      
      #res <- OLS_own(Y,X,0)$residuals
      #print(t(res)%*%res)
      #print(max(xsim))
      #n  <- length(Y)
      #k  <- ncol(X)
      #df <- n-k
      #sigma2 <- as.vector(t(res)%*%res)/df
      
      #With GLS I can use the assumption of alpha and sigma2
      diagonal <- sigma2*xsim^alpha 
      #I create a matrix of 0 of dimension of Res2
      
      omega_1 <- matrix(0,T,T) 
      #the diagonal is 1/sigma_n^2
      diag(omega_1) <- 1/diagonal
      
      GLS_static = GLS_own(Y,X,omega_1)
      
      
      beta_0_GLS[i] <- GLS_static[1,1]
      beta_1_GLS[i] <- GLS_static[2,1]
      
      stdvs_0[i] <- GLS_static[1,2]
      stdvs_1[i] <- GLS_static[2,2]
      
      ttest_matrix[i,j] <- (beta_1_GLS[i] - beta1_test[j])/stdvs_1[i] 
    }
    
#  }
  
  colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)
  
  
  beta_0_est <- mean(beta_0_GLS)
  beta_1_est <- mean(beta_1_GLS)
  
  #Numerical standard errors

  stdvs_0_num <- sd(beta_0_GLS)
  stdvs_1_num <- sd(beta_1_GLS)
  
  #analytical standard errors
  
  e <- rnorm(T,0,sigma2*xsim^alpha)
  #now I can have my y !
  Y <- X%*%beta + e

  #res <- OLS_own(Y,X,0)$residuals
  #print(t(res)%*%res)
  #print(max(xsim))
  #n  <- length(Y)
  #k  <- ncol(X)
  #df <- n-k
  #sigma2 <- as.vector(t(res)%*%res)/df
  
  #With GLS I can use the assumption of alpha and sigma2
  diagonal <- sigma2*xsim^alpha 
  #I create a matrix of 0 of dimension of T
  omega_1 <- matrix(0,T,T) 
  #the diagonal is 1/sigma_n^2
  diag(omega_1) <- 1/diagonal
  

  #accept the sigma2 of the function
  GLS_static = GLS_own(Y,X,omega_1)
  #
  stdvs_0_ana <- GLS_static[1,2]
  stdvs_1_ana <- GLS_static[2,2]
  #
  #Choose sigma2 =1
  #sigma2 =1
  #
  x <- X
  sigma2 <- as.vector(t(res)%*%res)/df
  
  var_01_ana   <- sigma2 * solve(t(x) %*% omega_1 %*%x)
  ###
  ###
  #var_0_ana <- var_01_ana[1,1]
  #var_1_ana <- var_01_ana[2,2]
  ###
  ####the diagonal elements are the std of Betas
  stdvs_0_ana <- sqrt(var_01_ana[1,1])
  stdvs_1_ana <- sqrt(var_01_ana[2,2])
  #
  #
  ##estimated standard errors
  stdvs_0_est <- mean(stdvs_0)
  stdvs_1_est <- mean(stdvs_1)
  
  
  #Row of estimation
  table_beta0 <- cbind(b_0,beta_0_est,stdvs_0_ana,stdvs_0_num,stdvs_0_est)
  colnames(table_beta0) <- c("population","estimated beta","analytical","numerical","estimated")
  

  table_beta1 <- cbind(b_1,beta_1_est,stdvs_1_ana,stdvs_1_num,stdvs_1_est)
  colnames(table_beta1) <- c("population","estimated beta","analytical","numerical","estimated")
  
  
  
  ###
  #C Compute size and power
  ###
  
  #loop over the t-tests and give me the Critical values for each one
  colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)
  
  #let's do a matrix of critical values, alpha is 5%
  #the t statistics follows a student t with 2 degrees of freedom 
  
  CV_beta1 <- qt(p=.05/2, df=T-2, lower.tail=FALSE)
  
  rej_function <- function(ttest,cv){
    if (abs(ttest) > cv){
      return(1)
    }else{
      return(0)
    }
  }
  
  
  #initialise matrix of rejection 
  rej_matrix <- matrix(0,nrow(ttest_matrix),ncol(ttest_matrix))
  #Lets store the 1 and 0 of rejection in a matrix
  for (j in 1:5){
    rej_matrix[,j]<-sapply(ttest_matrix[,j],rej_function, cv= CV_beta1)
  }
  #this is the size, the mean of column 1 for which beta = 1
  size_beta1 <- mean(rej_matrix[,1])
  #there is a size bias if the critical values are wrong !
  
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
  
#  for (j in 1:length(beta1_test)) {
    
    for (i in 1:repl) {
      #stochastic X: X = rnorm(T,0,sigma2)
      #let's get some errors, we define sigma 2 =1 earlier
      e <- rnorm(T,0,sigma2*xsim^alpha)
      #now I can have my y !
      Y <- X%*%beta + e
      
      
      #For EGLS, we assume that Var(miu_{i,t}) is sigma_i^2
      
      n  <- length(Y)
      k  <- ncol(X)
      
      OLS_out <- OLS_own(Y,X,0) 
      
      res <- OLS_out$residuals
      #I take absolute value in log because otherwise taking log(X) does not work
      #for teacher no need to take log(X)
      x <- X
      y <- log(res^2)
      
      ## Run OLS
      xy     <- t(x)%*%y #indeed I need some kind of X_i
      xxi    <- solve(t(x)%*%x)
      coefs  <- as.vector(xxi%*%xy)
      sigma_hat   <- as.vector(x%*%coefs)
      
      #for teacher sigma_hat is ok
      #But we think we need to get the initial sigma_hat to get the variance
      sigma_hat <- exp(sigma_hat)
      
      #now we create the matrix and put sigma_hat as the diagonal of that matrix
      omega_hat <- matrix(0,length(res), length(res)) 
      #I put back the diagonal element in the diagnonals
      #I take absolute value to have positive values 
      diag(omega_hat) <- 1/sigma_hat
      
      GLS_static = GLS_own (Y,X,omega_hat)
      
      
      beta_0_EGLS[i] <- GLS_static[1,1]
      beta_1_EGLS[i] <- GLS_static[2,1]
      
      stdvs_0[i] <- GLS_static[1,2]
      stdvs_1[i] <- GLS_static[2,2]
      
      ttest_matrix[i,j] <- (beta_1_EGLS[i] - beta1_test[j])/stdvs_1[i] 
    }
    
#  }
  
  colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)
  
  
  beta_0_est <- mean(beta_0_EGLS)
  beta_1_est <- mean(beta_1_EGLS)
  
  #Numerical standard errors
    
  var_0_num <- var(beta_0_EGLS)
  var_1_num <- var(beta_1_EGLS)
  
  stdvs_0_num <- sqrt(var_0_num)
  stdvs_1_num <- sqrt(var_1_num)

  #Analytical standard errors  
  sigma2 <- 1
  e <- rnorm(T,0,sigma2*xsim^alpha)
  #now I can have my y !
  Y <- X%*%beta + e
  
  OLS_out <- OLS_own(Y,X,0) 
  
  res <- OLS_out$residuals
  
  res2<- res%*%t(res)
  #the non-diagnonal elemets of res2 need to have 0
  #I take away the diagonals
  #for teacher no need to take log(X)
  x <- X

  y <- log(diag(res2))
  
  ## Run OLS
  xy     <- t(x)%*%y #indeed I need some kind of X_i
  xxi    <- solve(t(x)%*%x)
  coefs  <- as.vector(xxi%*%xy)
  sigma_hat   <- as.vector(x%*%coefs)
  #for teacher sigma_hat is ok
  sigma_hat <- exp(sigma_hat)
  
  hist(sigma_hat)
  
  #now we create the matrix and put sigma_hat as the diagonal of that matrix
  
  omega_hat <- matrix(0,length(res), length(res)) 
  #I put back the diagonal element in the diagnonals
  #I take absolute value to only have positive sigma_hat
  diag(omega_hat) <- 1/sigma_hat

  GLS_static = GLS_own (Y,X,omega_hat)
 
  stdvs_0_ana <- GLS_static[1,2]
  stdvs_1_ana <- GLS_static[2,2]
  stdvs_0_ana
  stdvs_1_ana
  
  #I need to do it by hand because the formula computes a different sigma2 
  #than the one we are supposed to assume
  #x <- X
  #
  #
  #var_01_MC   <- sigma2 * solve(t(x) %*% omega_hat %*%x)
  #
  #
  #var_0_MC <- var_01_MC[1,1]
  #var_1_MC <- var_01_MC[2,2]
  #
  ##the diagonal elements are the std of Betas
  #stdvs_0_ana <- sqrt(var_01_MC[1,1])
  #stdvs_1_ana <- sqrt(var_01_MC[2,2])
  
  #estimated standard errors
  stdvs_0_est <- mean(stdvs_0)
  stdvs_1_est <- mean(stdvs_1)
  
  
  #Row of estimation
  table_beta0 <- cbind(b_0,beta_0_est,stdvs_0_ana,stdvs_0_num,stdvs_0_est)
  colnames(table_beta0) <- c("population","estimated beta","analytical","numerical","estimated")
  
  #first row of OLS
  table_beta1 <- cbind(b_1,beta_1_est,stdvs_1_ana,stdvs_1_num,stdvs_1_est)
  colnames(table_beta1) <- c("population","estimated beta","analytical","numerical","estimated")
  
  
  
  ###
  #C Compute size and power
  ###
  
  #loop over the t-tests and give me the Critical values for each one
  colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)
  
  #let's do a matrix of critical values, alpha is 5%
  #the t statistics follows a student t with 2 degrees of freedom 
  
  CV_beta1 <- qt(p=.05/2, df=T-2, lower.tail=FALSE)
  
  rej_function <- function(ttest,cv){
    if (abs(ttest) > cv){
      return(1)
    }else{
      return(0)
    }
  }
  
  #initialise matrix of rejection 
  rej_matrix <- matrix(0,nrow(ttest_matrix),ncol(ttest_matrix))
  #Lets store the 1 and 0 of rejection in a matrix
  for (j in 1:5){
    rej_matrix[,j]<-sapply(ttest_matrix[,j],rej_function, cv= CV_beta1)
  }
  #this is the size, the mean of column 1 for which beta = 1
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
  
  #so there is a problem because the true std (analytical and numerical) should
  #coincide. In all cases, we correct for small (dividing by sqrt of T) so it really
  #should be the same.
  #example of beta_1
  
  #EGLS is consistent but biased
  #as N goes big -> the variance of EGLS is smaller than OLS
  
  beta_1_mat[,3:4]
  
  
  
  
  
 