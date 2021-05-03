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
setwd("C:/Users/qiszhang/Documents/programming/myrepo") #Windows

## Set output file directory: change to your own directory
#I want the output directory to be the same 
output = getwd() 
#so now everything is in one directory

## Install required packages (code is only executed if the packages are not yet installed)
if(!require(plm)){install.packages("plm")} ## panel data package
#install.packages("expm")
install.packages(readxl)


## Load required packages
library(plm) 
library(readxl)
#library(expm)

## Load your own functions
source("Case4_FE estimator.R")
source("Case4_GMMd estimator.R")

####################################################################
#Case 4.1. Panel data estimator                             ########
####################################################################
# Read from file
data  <- read_excel("Data_Baltagi.xlsx")

lg <- function(x)c(NA, x[1:(length(x)-1)])

data$ln.C_it_1 <- unlist(tapply(data$`ln C_it`, data$state, lg))

data  <- na.omit(data)

#We create the matrix D of dummies for each variable

y <- as.matrix(data[,10])

x <- as.matrix(data[,11:14])


idx <- sort(unique(data$state))
state_i <- matrix(0, nrow = nrow(data), ncol = length(idx))

for (j in 1:length(idx)) { 
  state_i[,j] <- as.integer(data$state == idx[j])
}

idx <- sort(unique(data$year))
year_i <- matrix(0, nrow = nrow(data), ncol = length(idx))

for (j in 1:length(idx)) { 
  year_i[,j] <- as.integer(data$year == idx[j])
}

#Let's write the demeaning matrix by individuals
#I need the identity matrix
NT <- nrow(data)
I   <- matrix(0, NT, NT)
diag(I) <- 1

P_i <- state_i %*% solve(t(state_i) %*% state_i) %*% t(state_i)

P_T <- year_i %*% solve(t(year_i) %*% year_i) %*% t(year_i)

y_demean <- y - P_i %*% y

x_demean <- x - P_i %*% x


y_demean_it <- y_demean - P_T %*% y_demean

x_demain_it <- x_demean - P_T %*% x_demean


n  <- length(y)
k  <- ncol(x)
df <- n-k

## Run OLS
x <- x
y <- y

xy     <- t(x)%*%y
xxi    <- solve(t(x)%*%x) #this is (X' X)^(-1)
coefs  <- as.vector(xxi%*%xy)

yhat   <- as.vector(x%*%coefs)
res    <- y-yhat
#sigma2 <- as.vector(t(res)%*%res/df)
sigma2 <- as.vector(t(res)%*%res)/df

#case where we only want to do 

  var <- sigma2*xxi
  stdvs <- sqrt(diag(var))
  
  tstats <- coefs/stdvs
  pvals  <- 2*(1-pt(abs(tstats),df))
  
  ## Save output
  names(coefs) <- colnames(x)
  
  coefs  <- round(coefs,3)
  stdvs  <- round(stdvs,3)
  tstats <- round(tstats,3)
  pvals  <- round(pvals,3)
  
  out = rbind(coefs,stdvs,tstats,pvals)
  out = t(out)

##First let's first difference our data!
  
data <- transform(data, dlnC_it = ave(`ln C_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnP_it = ave(`ln.P_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnPn_it = ave(`ln.Pn_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnY_it = ave(`ln.Y_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnC_it_1 = ave(ln.C_it, state, FUN = function(x) c(NA, diff(x))))

data<-na.omit(data)

#Extract T and N

N <-   length(unique(data$state))
T <-   length(unique(data$year))


y_fd <- as.matrix(data[,15])
x_fd <- as.matrix(data[,16:19])
  
#Let's create our H matrix

diagonal <- 2
offdiagonal<- -1
H <- matrix(0,T-2,T-2)
diag(H) <- diagonal
diag(H[-1,])<-offdiagonal
diag(H[,-1])<-offdiagonal
H

#Let's create the Zi matrix

yit <- seq(1,T-2)
n_inst <- sum(yit)
Z_i <- matrix(0,T-2,n_inst)
column <-1 

for (i in (1:(T-2))) {
  print(i)
  chunk <- as.numeric(yit[1:i])
  column <- column + i -1
  
  for (j in (1:length(chunk))) {
    Z_i[i,column+j-1] <- chunk[j]
    }
  
}



  
  