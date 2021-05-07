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
#it's fiiine !

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

colnames(data)[1] <- c('state')

lg <- function(x)c(NA, x[1:(length(x)-1)])

data$ln.C_it_1 <- unlist(tapply(data$`ln.C_it`, data$state, lg))

data  <- na.omit(data)

#We create the matrix D of dummies for each variable

y_S <- as.matrix(data[,10])

x_S <- as.matrix(data[,11:14])

#ok now we've got all the data we need!

####################################################################
###                Run OLS on static model                       ###
####################################################################

#My OLS function returns 2 things in a list

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
#the non-diagnonal elemets of res2 need to have 0
#I take away the diagonals
diagonal <- diag(res2) 
#I create a matrix of 0 of diam
omega_1 <- matrix(0,nrow(res2), ncol(res2)) 
#I put back the diagonal element in the diagnonal  
diag(omega_1) <- 1/diagonal


GLS_static = GLS_own (y_S,x_S,omega_1)


####################################################################
###                Run EGLS on static model                       ###
####################################################################

res <- OLS_static$residuals
k <- OLS_static$param

#We have 30 error terms per states
t <- 29

#There is one different sigma per state (46 states)
sigma_i <- seq(1:46)


#let's compute the sigma_i 
for (i in 1:(length(res)/t)){
  sigma_i[i] <- (t(res[(1+(i-1)*t):(i*t)])%*%(res[(1+(i-1)*t):(i*t)]))/(t-k) 
}
#Maybe K is 46 but not sure

#Now we create the diagonal of sigma_hat

sigma_est<-rep(sigma_i,each=t)

#now we create the matrix and put sigma_hat as the diagonal of that matrix

omega_hat <- matrix(0,length(res), length(res)) 
#I put back the diagonal element in the diagnonals
diag(omega_hat) <- sigma_est

omega_hat

GLS_static = GLS_own (y_S,x_S,omega_hat)



##################################################################
###     Run OLS on dynamic model in first difference           ###
##################################################################
data    = read.table("Data_Baltagi.csv", header = TRUE,sep = ";")

colnames(data)[1] <- c('state')

lg <- function(x)c(NA, x[1:(length(x)-1)])

data$ln.C_it_1 <- unlist(tapply(data$`ln.C_it`, data$state, lg))

data  <- na.omit(data)

#We create the matrix D of dummies for each variable

y <- as.matrix(data[,10])

x <- as.matrix(data[,11:14])

##First let's first difference our data!

data <- transform(data, dlnC_it = ave(`ln.C_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnP_it = ave(`ln.P_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnPn_it = ave(`ln.Pn_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnY_it = ave(`ln.Y_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnC_it_1 = ave(ln.C_it_1, state, FUN = function(x) c(NA, diff(x))))

data<-na.omit(data)


#Extract T and N

N <-   length(unique(data$state))

#we loose an obs with lag and another one with difference

y_fd <- as.matrix(data[,15])
x_fd <- as.matrix(data[,16:19])


OLS_dynamic = OLS_own(y_fd,x_fd,1)

##################################################################
###     Run EGLS on dynamic model in first difference           ###
##################################################################



OLS_dynamic = OLS_own(y_fd,x_fd,0)

res <- OLS_dynamic$residuals

k <- OLS_dynamic$param

#We have 29 error terms per states because we do first difference
t <- 28

#There is one different sigma per state (46 states)
sigma_i <- seq(1:46)


#let's compute the sigma_i 
for (i in 1:(length(res)/t)){
  sigma_i[i] <- (t(res[(1+(i-1)*t):(i*t)])%*%(res[(1+(i-1)*t):(i*t)]))/(t-k) 
}
#Maybe K is 46 but not sure

#Now we create the diagonal of sigma_hat

sigma_est<-rep(sigma_i,each=t)

#now we create the matrix and put sigma_hat as the diagonal of that matrix

omega_hat <- matrix(0,length(res), length(res)) 
#I put back the diagonal element in the diagnonals
diag(omega_hat) <- sigma_est

omega_hat


OLS_dynamic = GLS_own(y_D_diff,x_D_diff,sigma_est)

