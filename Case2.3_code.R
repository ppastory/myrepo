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
repl <- 1000 #less number of replication to work on the code
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

brepl <- 50

#I am initialising the vectors in which beta and other stuff will arrive
beta_0_0LSbp <- rep(0,brepl)
beta_1_OLSbp <- rep(0,brepl)

#I am initialising the vectors in which beta and other stuff will arrive
beta_0_0LSbw <- rep(0,brepl)
beta_1_OLSbw <- rep(0,brepl)



#grid of beta1 
beta1_test <- as.matrix(t(c(1,0.95,0.90,0.75,0.5)))

#ttest matrix with column are different beta1 and 
#the rows are one ttest per simulation

ttest_matrix <- matrix(0,repl, length(beta1_test))


data <- cbind(Y,X)


for (j in 1:length(beta1_test)) {
  
  for (i in 1:brepl) {

    boot_pair = data[unlist(sample(as.data.frame(matrix(1:nrow(data),nrow = 2)),100,replace=T)),]
    
    Ybp<-boot_pair[,1]
    Xbp<-boot_pair[,-1]
    
    OLS_out <- OLS_own(Ybp,Xbp,0) 
    
    
    beta_0_0LSbp[i] <- OLS_out$estimation[1,1]
    beta_1_OLSbp[i] <- OLS_out$estimation[2,1]
    
    #stdvs_0[i] <- OLS_out$estimation[1,2]/sqrt(T)
    #stdvs_1[i] <- OLS_out$estimation[2,2]/sqrt(T)
    #
    #ttest_matrix[i,j] <- (beta_1_OLS[i] - beta1_test[j])/stdvs_1[i] #divided by sqrt T
    
    shock <- sample(c(-1,1), replace=TRUE, size=200)
    errors_b <- shock*e
    
    #now I can have my y !
    Y_b <- X%*%beta + errors_b
    
    OLS_out <- OLS_own(Y_b,X,0) 
    
    
    beta_0_0LSbw[i] <- OLS_out$estimation[1,1]
    beta_1_OLSbw[i] <- OLS_out$estimation[2,1]
    
  }
  
}

colnames(ttest_matrix) <- c(1,0.95,0.90,0.75,0.5)

beta_0_barp <- mean(beta_0_0LSbp)
beta_1_barp <- mean(beta_1_OLSbp)
#let's get the numerical standard errors
#let's get the numerical standard errors -> truuuue

var_0_nump <- var(beta_0_0LSbp)/T ##divide by T
var_1_nump <- var(beta_1_OLSbp)/T

stdvs_0_nump <- sqrt(var_0_nump)
stdvs_1_nump <- sqrt(var_1_nump)

beta_0_barw <- mean(beta_0_0LSbw)
beta_1_barw <- mean(beta_1_OLSbw)
#let's get the numerical standard errors
#let's get the numerical standard errors -> truuuue

var_0_numw <- var(beta_0_0LSbw)/T ##divide by T
var_1_numw <- var(beta_1_OLSbw)/T

stdvs_0_numw <- sqrt(var_0_numw)
stdvs_1_numw <- sqrt(var_1_numw)

#first row of OLS
table_std_num <- cbind(stdvs_0_nump,stdvs_1_nump,stdvs_0_numw,stdvs_1_numw)
colnames(table_std_num) <- c("pair std beta0","pair std beta1","wild std beta0","wild std beta1")
table_std_num






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
