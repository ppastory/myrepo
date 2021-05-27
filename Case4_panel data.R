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
data    = read.table("Data_Baltagi.csv", header = TRUE,sep = ";")

colnames(data)[1] <- c('state')

lg <- function(x)c(NA, x[1:(length(x)-1)])

data$ln.C_it_1 <- unlist(tapply(data$`ln.C_it`, data$state, lg))

data  <- na.omit(data)

#We create the matrix D of dummies for each variable

FE_est = FE(data,1)




###################################################
############Restrict the lag length###############
##################################################


#let's say that I want max 3 lags
max_lag <- 5

el <- max_lag-1
#Let's loop over the states and create our big matrix

n_inst <- sum(seq(1,(max_lag-1))) + ((T-2)-length(seq(1,(max_lag-1))))*(max_lag-1)


#The chunks have size 27 !!!
#because yi27 is the last instrument of delta_ei29
y_1 <- y[1:(T-2)]

Z_1 <- matrix(0,T-2,n_inst)
column <-1 

for (i in (1:(T-2))) {
  
  #Chunk is the bunch of yi that we are going to put in the Z_i matrix
  chunk <- as.numeric(y_1[1:i])
  #print(chunk)
  
  if (length(chunk) < (max_lag-1)){
    for (j in 1:length(chunk)){
      Z_1[i,column+j-1] <- chunk[j]
    }
    # Z_1[i,column] <- chunk
    column <- column + i
    
  } else {
    
    print(column)
    el <- max_lag-1
    
    for (j in 1:(max_lag-1)){
      Z_1[i,column+j-1] <- tail(chunk,el)[j]
    }
    
    column <- column + (max_lag-1)
  }
}


Z <- Z_1

#We start at the column 30 and go by chunks of 29
#but we will only loop over the Y chunks until 27 !
#that is because yi27 will instrument delta_ei29
for (sst in seq(30, length(y), T)) {
  
  print(sst)
  #The second chunk for example goes from 30 to 
  #30 + 27 
  y_i <- y[sst:(sst+27)]
  
  #The matrix of instruments has size 27 x 378 ! 
  #each time period we use more lags
  Z_i <- matrix(0,T-2,n_inst)
  column <-1 
  
  for (i in (1:(T-2))) {
    
    #Chunk is the bunch of yi that we are going to put in the Z_i matrix
    chunk <- as.numeric(y_1[1:i])
    #print(chunk)
    
    if (length(chunk) < (max_lag-1)){
      for (j in 1:length(chunk)){
        Z_i[i,column+j-1] <- chunk[j]
      }
      # Z_1[i,column] <- chunk
      column <- column + i
      
    } else {
      
      #print(column)
      el <- max_lag-1
      
      for (j in 1:(max_lag-1)){
        Z_i[i,column+j-1] <- tail(chunk,el)[j]
      }
      
      column <- column + (max_lag-1)
    }
  }
  
  
  Z <- rbind(Z,Z_i)
  
}

#changing number of instrumeeents!
n_inst_var <- N*(T-2)

#Imagine our big H is very big -> we try to create 
diagonal <- 2
offdiagonal<- -1
H <- matrix(0,n_inst_var,n_inst_var)
diag(H) <- diagonal
diag(H[-1,])<-offdiagonal
diag(H[,-1])<-offdiagonal
H

#Now we need to create the big X and Y matrix but it is not simply x and y

#Remember we can only start instrumenting in period 3 until period 29 
#this means that we take the 27 last periods for every state!
#for every state I need to kill the first 2 periods
data_s <- data

data_s[data_s == 65] <- NA
data_s <- na.omit(data_s)

y_fd_s <- as.matrix(data_s[,15])
x_fd_s <- as.matrix(cbind(Cnst=1,data_s[,16:19]))

#we have all the ingredients we need to compute our gamma
big_Z <- as.matrix(cbind(cbind(Cnst=1,x_fd_s[,1],x_fd_s[,2],x_fd_s[,3],Z)))

Z <- big_Z

W_notinv <- t(Z) %*% H %*% Z 
#we have the big W optimal
W_opt <- solve(W_notinv, tol = 1e-20)


gamma <- solve(t(x_fd_s) %*% Z %*% W_opt %*% t(Z) %*% x_fd_s) %*%  t(x_fd_s) %*% Z %*% W_opt %*% t(Z) %*% y_fd_s

yhat   <- as.vector(x_fd_s%*%gamma)
res    <- y-yhat

sigma2 <- as.vector(t(res)%*%res)/2

var_gamma <- sigma2* solve(t(x_fd_s) %*% Z %*% W_opt %*% t(Z) %*% x_fd_s)

var_g <- diag(var_gamma)

std_P <- sqrt(var_g[1])
std_Pn <- sqrt(var_g[2])
std_Y <- sqrt(var_g[3])
std_Ct_1 <- sqrt(var_g[4])
