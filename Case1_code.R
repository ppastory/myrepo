## Sample R script 
## Using the pooled OLS estimator to estimate the Cigarette demand example from the Baltagi textbook

#version 1.4

####################################################################
###              Set-up for use                                  ###
####################################################################
rm(list = ls())   # Clear workspace

## Set input file directory

#let's check what is the current directory
getwd()

#We don't have to sent the current directory because it is our directory myrepo
#C:/Users/ppastory/Documents/programming/myrepo
setwd("C:\\Users\\qiszhang\\Documents\\programming/myrepo") #Windows

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
###########        Load data into the workspace          ###########
####################################################################
# Read from file
data    = read.table("Data_Baltagi.csv", header = TRUE,sep = ";")

#rename column
colnames(data)[1] <- c('state')

#function to take the lag
lg <- function(x)c(NA, x[1:(length(x)-1)])

#take the lag
data$ln.C_it_1 <- unlist(tapply(data$`ln.C_it`, data$state, lg))

#take out the na 
data  <- na.omit(data)

#Create the data set static

y_S <- as.matrix(data[,10])

x_S <- as.matrix(cbind(Cnst=1, data[,11:14]))


####################################################################
###                Run OLS on static model                       ###
####################################################################

#My OLS function returns the estimation, the residuals and the parameters

OLS_static = OLS_own(y_S,x_S,0)

####################################################################
###         Run OLS on static model with White correction        ###
####################################################################

OLS_static_white_errors = OLS_own(y_S,x_S,1)


####################################################################
###                Run GLS on static model                       ###
####################################################################

res <- OLS_static$residuals

### I compute Omega outside of the function
res2       <- res%*%t(res)
#the non-diagonal elements of res2 need to have 0
#I take away the diagonals
diagonal <- diag(res2) 
#I create a matrix of 0 of diag
P <- matrix(0,nrow(res2), ncol(res2)) 

diag(P) <- sqrt(diagonal)

omega_1bis <- t(P) %*% P

GLS_static <- GLS_own(y_S,x_S,omega_1bis)

####################################################################
###                Run EGLS on static model                       ###
####################################################################

res <- OLS_static$residuals
k <- OLS_static$param

#We have 29 error terms per states
t <- 29

#There is one different sigma per state (46 states)
sigma2_i <- seq(1:46)


#let's compute the sigma_i 
for (i in 1:(length(res)/t)){
  sigma2_i[i] <- (t(res[(1+(i-1)*t):(i*t)])%*%(res[(1+(i-1)*t):(i*t)]))/(t-k) 
}


#Now we create the diagonal of sigma_hat
sigma2_est<-rep(sigma2_i,each=t)

#now we create the matrix and put sigma_hat as the diagonal of that matrix
omega_hat <- matrix(0,length(res), length(res)) 

#I put back the diagonal element in the diagnonals
diag(omega_hat) <- sigma2_est

#run GLS
GLS_static = GLS_own (y_S,x_S,omega_hat)


##################################################################
###     Run OLS on dynamic model in first difference           ###
##################################################################
data    = read.table("Data_Baltagi.csv", header = TRUE,sep = ";")

colnames(data)[1] <- c('state')

lg <- function(x)c(NA, x[1:(length(x)-1)])

data$ln.C_it_1 <- unlist(tapply(data$`ln.C_it`, data$state, lg))

data  <- na.omit(data)

#extract the data we need

y <- as.matrix(data[,10])

x <- as.matrix(data[,11:14])

##First let's first difference our data!

data <-   transform(data, dlnC_it = ave(`ln.C_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnP_it = ave(`ln.P_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnPn_it = ave(`ln.Pn_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnY_it = ave(`ln.Y_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnC_it_1 = ave(ln.C_it_1, state, FUN = function(x) c(NA, diff(x))))

data<-na.omit(data)


#Extract T and N

N <-   length(unique(data$state))

#we loose an obs with lag and another one with difference pre state

y_fd <- as.matrix(data[,15])
x_fd <- as.matrix(cbind(Cnst=1,data[,16:19]))


OLS_dynamic = OLS_own(y_fd,x_fd,0)

##################################################################
###     Run EGLS on dynamic model in first difference           ###
##################################################################

OLS_dynamic = OLS_own(y_fd,x_fd,0)
#extract errors
res <- OLS_dynamic$residuals

#extract parameters
k <- OLS_dynamic$param

#We have 29 error terms per states because we do first difference, there are 28
t <- 28

#There is one different sigma per state (46 states)
sigma_i <- seq(1:46)


#let's compute the sigma_i 
for (i in 1:(length(res)/t)){
  sigma_i[i] <- (t(res[(1+(i-1)*t):(i*t)])%*%(res[(1+(i-1)*t):(i*t)]))/(t-k) 
}


#Now we create the diagonal of sigma_hat

sigma_est<-rep(sigma_i,each=t)

#now we create the matrix and put sigma_hat as the diagonal of that matrix

omega_hat <- matrix(0,length(res), length(res)) 
#I put back the diagonal element in the diagnonals
diag(omega_hat) <- sigma_est

omega_hat


OLS_dynamic = GLS_own(y_fd,x_fd,omega_hat)

