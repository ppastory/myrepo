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
lnC_it  = data$ln.C_it
lnP_it  = data$ln.P_it
lnPn_it = data$ln.Pn_it
lnY_it  = data$ln.Y_it

# Construct matrices (for static model)
y_S<-as.vector(lnC_it)
x_S<-as.matrix(cbind(Cnst=1,lnP_it,lnPn_it,lnY_it))

# Construct matrices (for dynamic model)
pdt <- pdata.frame(data)            # use plm package to construct lag of lnC_it
pdt$ln.C_itL1 <- lag(pdt$ln.C_it)
pdt_noNaN     <- na.omit(pdt)

lnC_it_D    <- pdt_noNaN$ln.C_it
lnC_itL1_D  <- pdt_noNaN$ln.C_itL1
lnP_it_D    <- pdt_noNaN$ln.P_it
lnPn_it_D   <- pdt_noNaN$ln.Pn_it
lnY_it_D    <- pdt_noNaN$ln.Y_it

y_D<-as.vector(lnC_it_D)
x_D<-as.matrix(cbind(Cnst=1,lnC_itL1_D,lnP_it_D,lnPn_it_D,lnY_it_D))

# Take first differences
lnC_it_D_diff   = diff(lnC_it_D, differences = 1)
lnC_itL1_D_diff = diff(lnC_itL1_D, differences = 1)
lnP_it_D_diff   = diff(lnP_it_D, differences = 1)                                    # Take first differences                
lnPn_it_D_diff  = diff(lnPn_it_D, differences = 1)
lnY_it_D_diff   = diff(lnY_it_D, differences = 1)

y_D_diff<-as.vector(lnC_it_D_diff)
x_D_diff<-as.matrix(cbind(Cnst=1,lnC_itL1_D_diff,lnP_it_D_diff,lnPn_it_D_diff,lnY_it_D_diff))

y_D_diff = na.omit(y_D_diff) #Remove NA's otherwise the OLS functions returns an error
x_D_diff = na.omit(x_D_diff)

#ok now we've got all the data we need!

####################################################################
###                Run OLS on static model                       ###
####################################################################

#My OLS function returns 2 things in a list

OLS_static = OLS_own(y_S,x_S,0)
OLS_static$residuals

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
res2 <- matrix(0,nrow(res2), ncol(res2)) 
#I put back the diagonal element in the diagnonals
#but this time I take the inverse to get P
diag(res2) <- 1/sqrt(diagonal)
#res2 is sigma2 omega in the formulas
P <- res2

omega_1    <- t(P)%*%P


GLS_static = GLS_own (y_S,x_S,omega_1)


####################################################################
###                Run EGLS on static model                       ###
####################################################################

res <- OLS_static$residuals
k <- OLS_static$param

#We have 30 error terms per states
t <- 30

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


OLS_dynamic = OLS_own(y_D_diff,x_D_diff,1)
OLS_dynamic = GLS_own(y_D_diff,x_D_diff,sigma_est)

