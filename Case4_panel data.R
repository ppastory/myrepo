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

  T <-   length(unique(data$year)) #the real T
  
  
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
  
#Let's create the Zi matrix
#we need to start with the Z_1 matrix that below which we are going to stack the other Z_i


#Let's loop over the states and create our big matrix

yit <- seq(1,T-2)
n_inst <- sum(yit)

#The chunks have size 27 !!!
#because yi27 is the last instrument of delta_ei29
y_1 <- y[1:(T-2)]

#I need to take the X_fd_s, the X that loose the first 2 observation by individuals
#because you do FD and because you need a lag (we start in delta xi3 basically)

data_s <- data

data_s[data_s == 65] <- NA
data_s <- na.omit(data_s)

y_fd_s <- as.matrix(data_s[,15])

#we what the order of x_fd_s to be the same as in the Z_i matrix
#and as in the book
x_fd_s <- as.matrix(data_s[,16:19])

x_fd_s[,1] <- as.matrix(data_s[,19])
x_fd_s[,2] <- as.matrix(data_s[,16])
x_fd_s[,3] <- as.matrix(data_s[,18])
x_fd_s[,4] <- as.matrix(data_s[,17])
#We check the names and they are wrong so let's change the column names
colnames(x_fd_s)[1] <- colnames(x_fd_s)[4]
colnames(x_fd_s)[2] <- c('dlnP_it')
colnames(x_fd_s)[3] #it's ok
colnames(x_fd_s)[4] <- c('dlnPn_it')
x_fd_s

#The number of columns of Z_1 is we would have without adding is 1 + 2 + .. + 27 = ninst
#but now for each chunk 1 + 2 + ... + 27, we had the 3 explanatory variables as instrument

Z_1 <- matrix(0,T-2,n_inst)
column <-1 

for (i in (1:(T-2))) {
  
  #Chunk is the bunch of yi that we are going to put in the Z_i matrix
  chunk <- as.numeric(y_1[1:i])
  #print(chunk)
  column <- column + i -1
  
  for (j in (1:length(chunk))) {
    Z_1[i,column+j-1] <- chunk[j]
  }
  
}


#Z_1 

Z <- Z_1

#We start at the column 30 and go by chunks of 29
#but we will only loop over the Y chunks until 27 !
#that is because yi27 will instrument delta_ei29
for (sst in seq(30, length(y), T)) {

  #print(sst)
  #The second chunk for example goes from 30 to 
  #30 + 27 
  y_i <- y[sst:(sst+27)]
  
  #The matrix of instruments has size 27 x 378 ! 
  #each time period we use more lags
  Z_i <- matrix(0,T-2,n_inst)
  column <-1 

  for (i in (1:(T-2))) {
    #print(i)
    chunk <- as.numeric(y_i[1:i])
    column <- column + i -1
  
    for (j in (1:length(chunk))) {
      Z_i[i,column+j-1] <- chunk[j]
    }
  
  }


Z <- rbind(Z,Z_i)

}

n_inst_var <- N*(T-2)

#Imagine our big H is very big -> we try to create 
diagonal <- 2
offdiagonal<- -1
H <- matrix(0,n_inst_var,n_inst_var)
diag(H) <- diagonal
diag(H[-1,])<-offdiagonal
diag(H[,-1])<-offdiagonal
H
#this is a matrix with diagonal 2 and -1 on each side 


W_notinv <- t(Z) %*% H %*% Z 
#we have the big W optimal
W_opt <- solve(W_notinv)

#Now we need to create the big X and Y matrix but it is not simply x and y

#Remember we can only start instrumenting in period 3 until period 29 
#this means that we take the 27 last periods for every state!
#for every state I need to kill the first 2 periods
data_s <- data

data_s[data_s == 65] <- NA
data_s <- na.omit(data_s)

y_fd_s <- as.matrix(data_s[,15])
x_fd_s <- as.matrix(data_s[,16:19])

#we have all the ingredients we need to compute our gamma
big_Z <- as.matrix(cbind(x_fd_s[,1],x_fd_s[,2],x_fd_s[,3],Z))

Z <- big_Z

W_notinv <- t(Z) %*% H %*% Z 
#we have the big W optimal
W_opt <- solve(W_notinv, tol = 1e-20)

gamma <- solve(t(x_fd_s) %*% Z %*% W_opt %*% t(Z) %*% x_fd_s) %*%  t(x_fd_s) %*% Z %*% W_opt %*% t(Z) %*% y_fd_s

y <- y_fd_s
yhat   <- as.vector(x_fd_s%*%gamma)
res    <- y-yhat

sigma2 <- as.vector(t(res)%*%res)/2

var_gamma <- sigma2* solve(t(x_fd_s) %*% Z %*% W_opt %*% t(Z) %*% x_fd_s)

var_g <- diag(var_gamma)

std_P <- sqrt(var_g[1])
std_Pn <- sqrt(var_g[2])
std_Y <- sqrt(var_g[3])
std_Ct_1 <- sqrt(var_g[4])

###################################################
############Restrict the lag length###############
##################################################


#let's say that I want max 3 lags
max_lag <- 4

el <- max_lag-1
#Let's loop over the states and create our big matrix

n_inst <- max_lag*(T-2)

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
    
    for (j in 1:length(chunk)){
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
    chunk <- as.numeric(y_i[1:i])
    #print(chunk)
    
    if (length(chunk) < 2){
      Z_i[i,column] <- chunk
      column <- column + 1
      
    } else {
      
      el <- max_lag-1
      
      Z_i[i,column] <- tail(chunk,el)[1]
      Z_i[i,column+1] <- tail(chunk,el)[2]
      column <- column + 2
    }
    
  }
  
  
  Z <- rbind(Z,Z_i)
  
}

n_inst_var <- N*(T-2)

#Imagine our big H is very big -> we try to create 
diagonal <- 2
offdiagonal<- -1
H <- matrix(0,n_inst_var,n_inst_var)
diag(H) <- diagonal
diag(H[-1,])<-offdiagonal
diag(H[,-1])<-offdiagonal
H


W_notinv <- t(Z) %*% H %*% Z 
#we have the big W optimal
W_opt <- solve(W_notinv)

#Now we need to create the big X and Y matrix but it is not simply x and y

#Remember we can only start instrumenting in period 3 until period 29 
#this means that we take the 27 last periods for every state!
#for every state I need to kill the first 2 periods
data_s <- data

data_s[data_s == 65] <- NA
data_s <- na.omit(data_s)

y_fd_s <- as.matrix(data_s[,15])
x_fd_s <- as.matrix(data_s[,16:19])

#we have all the ingredients we need to compute our gamma
big_Z <- as.matrix(cbind(x_fd_s[,1],x_fd_s[,2],x_fd_s[,3],Z))

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
