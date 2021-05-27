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
#it's fine

#We don't have to sent the current directory because it is our directory myrepo
#C:/Users/ppastory/Documents/programming/myrepo
setwd("C:/Users/ppastory/Documents/programming/myrepo") #Windows

## Set output file directory: change to your own directory
#I want the output directory to be the same 
output = getwd() 
#so now everything is in one directory

## Install required packages (code is only executed if the packages are not yet installed)
#if(!require(plm)){install.packages("plm")} ## panel data package
#install.packages("xslx")
#install.packages(readxl)
#install.packages("dplyr")                       # Install dplyr package

## Load required packages
library(plm) 
library(readxl)
library(xslx)
library(dplyr)                                # Load dplyr package


## Load your own functions
source("Case4_FE estimator.R")
source("Case4_panel_GMMd_function.R")

####################################################################
#Case 4.2. System GMM (first style)                         ########
####################################################################
# Read from file
data    = read.table("Data_Baltagi.csv", header = TRUE,sep = ";")

colnames(data)[1] <- c('state')

lg <- function(x)c(NA, x[1:(length(x)-1)])

data$ln.C_it_1 <- unlist(tapply(data$`ln.C_it`, data$state, lg))


#we extract the data we need

y <- as.matrix(data[,10])

x <- as.matrix(data[,11:14])



##First let's first difference our data!

data <-   transform(data, dlnC_it = ave(`ln.C_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnP_it = ave(`ln.P_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnPn_it = ave(`ln.Pn_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnY_it = ave(`ln.Y_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnC_it_1 = ave(ln.C_it_1, state, FUN = function(x) c(NA, diff(x))))

#all the data has time fixed effect

var <- data[,c(2,15:19)]

#var %>% group_by(year) %>% mutate(Mean <- mean(var$dlnP_it, na.rm=TRUE))

var <- var %>%
  group_by(year)  %>%
  mutate(across(c("dlnC_it","dlnP_it",  "dlnPn_it","dlnY_it","dlnC_it_1"), ~ .x - mean(.x, na.rm=TRUE), .names = "tdm_{col}"))

#replace value with the demean value 

data[,15:19] <- var[,7:11]



##Full instrument matrix
output <- dGMM(data,1,0,0)
output

##Stacked matrix
output <- dGMM(data,1,1,0)
output

##Max lag matrix
#put min max lag 2 because the minimum lag is 2
output <- dGMM(data,0,0,10)
output

#lag 1 has a problem 

#Extract T and N

N <-   length(unique(data$state))
T <-   length(unique(data$year)) #the real T


#Let's loop over the states and create our big matrix

yit <- seq(1,T-2)
n_inst <- sum(yit)

#The chunks have size 28 !!!
#because yi27 is the last instrument of delta_ei29
y_1 <- data[1:(T-2),10]



#The number of columns of Z_1 is we would have without adding is 1 + 2 + .. + 27 = ninst
#but now for each chunk 1 + 2 + ... + 27, we had the 3 explanatory variables as instrument

Z_1 <- matrix(0,T-2,n_inst+((T-2)*3))
column <-2 

for (i in (1:(T-2))) {
  
  #Chunk is the bunch of yi that we are going to put in the Z_i matrix
  chunk <- as.numeric(y_1[1:i])
  #print(chunk)
  exog_var <- data[i+2,16:18]
  
  endog_exog <- append(as.numeric(y_1[1:i]),as.numeric(exog_var))

    
  for (j in (1:length(endog_exog))) {
    Z_1[i,(column+j-2)] <- endog_exog[j]
  }
  column <- column + i +3
}


x_1 <- exog_var

x_i <- x_1


Z_i <- matrix(0,T-2,n_inst+((T-2)*3))

Z <- Z_1


#We start at the column 30 and go by chunks of 29
#but we will only loop over the Y chunks until 27 !
#that is because yi27 will instrument delta_ei29
for (sst in seq(31, nrow(data), T)) {
  
  #print(sst)
  #The second chunk for example goes from 30 to 
  #30 + 27 
  y_i <- data[sst:(sst+27),10]
  
  x_i <- data[(sst+2):((sst+2)+28),16:18]
  
  #The matrix of instruments has size 27 x 378 ! 
  #each time period we use more lags
  Z_i <- matrix(0,T-2,n_inst+((T-2)*3))
  
  column <-2 
  #Next time I will start in row 57 which is row 59 -2 
  #and so one the keeps increase by 1 for each loop, this -2 takes the name iter
  
  
  for (i in (1:(T-2))) {
    
    #Chunk is the bunch of yi that we are going to put in the Z_i matrix
    chunk <- as.numeric(y_1[1:i])
    #print(chunk)
    exog_var <- x_i[i,]
    
    
    endog_exog <- append(as.numeric(y_1[1:i]),as.numeric(exog_var))

    for (j in (1:length(endog_exog))) {
      Z_i[i,(column+j-2)] <- endog_exog[j]
    }
    column <- column + i +3
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

#this is a matrix with diagonal 2 and -1 on each side 

data_reg <- na.omit(data)

x1 <- as.matrix(na.omit(data_reg[,16:18]))
x2 <- as.matrix(na.omit(data_reg[,19]))
x<- cbind(x2,x1) #V1 like vendogemous

data_reg <- na.omit(data)
y <- as.matrix(na.omit(data_reg[,15]))

W_notinv <- t(Z) %*% H %*% Z 

#we have the big W optimal
W_opt <- solve(W_notinv,tol = 1e-22)


gamma <- solve(t(x) %*% Z %*% W_opt %*% t(Z) %*% x) %*%  t(x) %*% Z %*% W_opt %*% t(Z) %*% y

y_pred <- y
yhat   <- as.vector(x%*%gamma)
res    <- y_pred-yhat

n  <- length(y)
k  <- ncol(x)
df <- n-k

sigma2 <- as.vector(t(res)%*%res)/df

sigma2 <- sigma2/2

var_gamma <- sigma2* solve(t(x) %*% Z %*% W_opt %*% t(Z) %*% x)


#std_P <- sqrt(var_g[1])
#std_Pn <- sqrt(var_g[2])
#std_Y <- sqrt(var_g[3])
#std_Ct_1 <- sqrt(var_g[4])


stdvs_BGMM_sys  <-  sqrt(diag(var_gamma))
stdvs_BGMM_sys

tstats_BGMM_sys <- gamma/stdvs_BGMM_sys
tstats_BGMM_sys <- c(tstats_BGMM_sys)

n  <- length(y)
k  <- ncol(x)
df <- n-k

pvals_GGMM_sys   <- 2*(1-pt(abs(tstats_BGMM_sys),df))


#save output
names(gamma) <- colnames(x)
names(tstats_BGMM_sys) <- colnames(x)

#creating the table
coefs  <- round(gamma,3)
stdvs  <- round(stdvs_BGMM_sys,3)
tstats <- round(tstats_BGMM_sys,3)
pvals  <- round(pvals_GGMM_sys,3)


out_pvals_BGMM_sys = cbind(coefs,stdvs,tstats,pvals)


###################################################
############Restrict the lag length###############
##################################################


#let's say that I want max 5 lags
max_lag <- 10

#el is the number of maximum instruments per period
el <- max_lag-1
#Let's loop over the states and create our big matrix

n_inst <- sum(seq(1,(max_lag-1))) + ((T-2)-length(seq(1,(max_lag-1))))*(max_lag-1)


#The chunks have size 28 !!!
#because yi28 is the last instrument of delta_ei30

y_1 <- data[1:(T-2),10]

#number of instruments and exog variables entering the z matrix
n_var <- n_inst + (3*(T-2))

Z_1 <- matrix(0,T-2,n_var)
column <-2 

for (i in (1:(T-2))) {
  
  #Chunk is the bunch of yi that we are going to put in the Z_i matrix
  chunk <- as.numeric(y_1[1:i])
  #print(chunk)
  exog_var <- data[i+2,16:18]
  
  if (length(chunk) < (max_lag)){
    
    endog_exog <- append(as.numeric(y_1[1:i]),as.numeric(exog_var))
    
    for (j in 1:length(endog_exog)){
      print(endog_exog[j])
      Z_1[i,(column+j-2)] <- endog_exog[j]
      
    }
    
    column <- column + i +3
    
  } else {
    
    el <- max_lag-1
    
    endog_exog <- append(as.numeric(tail(chunk,el)),as.numeric(exog_var))
    
    
    for (j in 1:length(endog_exog)){
      
      Z_1[i,column+j-2] <- endog_exog[j]
    }
    
    column <- column + el + 3
  }
}


Z <- Z_1

#We start at the column 30 and go by chunks of 29
#but we will only loop over the Y chunks until 27 !
#that is because yi27 will instrument delta_ei29
for (sst in seq(31, nrow(data), T)) {
  
  #print(sst)
  #The second chunk for example goes from 30 to 
  #30 + 27 
  y_i <- data[sst:(sst+27),10]
  
  x_i <- data[(sst+2):((sst+2)+28),16:18]
  
  #The matrix of instruments has size 27 x 378 ! 
  #each time period we use more lags
  Z_i <- matrix(0,T-2,n_inst+((T-2)*3))
  
  column <-2 
  #Next time I will start in row 57 which is row 59 -2 
  #and so one the keeps increase by 1 for each loop, this -2 takes the name iter
  
  
  for (i in (1:(T-2))) {
    
    
    #Chunk is the bunch of yi that we are going to put in the Z_i matrix
    chunk <- as.numeric(y_i[1:i])
    #print(chunk)
    exog_var <- x_i[i,]
    
    
    if (length(chunk) < (max_lag)){
      
      endog_exog <- append(as.numeric(y_i[1:i]),as.numeric(exog_var))
      
      for (j in 1:length(endog_exog)){
        
        Z_i[i,(column+j-2)] <- endog_exog[j]
        
      }
      
      column <- column + i +3
      
    } else {
      
      el <- max_lag-1
      
      endog_exog <- append(as.numeric(tail(chunk,el)),as.numeric(exog_var))
      
      
      for (j in 1:length(endog_exog)){
        
        Z_i[i,column+j-2] <- endog_exog[j]
      }
      
      column <- column + el + 3
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



data_reg <- na.omit(data)

x1 <- as.matrix(na.omit(data_reg[,16:18]))
x2 <- as.matrix(na.omit(data_reg[,19]))
x<- cbind(x2,x1) #V1 like vendogemous

data_reg <- na.omit(data)
y <- as.matrix(na.omit(data_reg[,15]))

W_notinv <- t(Z) %*% H %*% Z 



#we have the big W optimal
W_opt <- solve(W_notinv,tol=1e-21)


gamma <- solve(t(x) %*% Z %*% W_opt %*% t(Z) %*% x) %*%  t(x) %*% Z %*% W_opt %*% t(Z) %*% y

y_pred <- y
yhat   <- as.vector(x%*%gamma)
res    <- y_pred-yhat

sigma2 <- as.vector(t(res)%*%res)/2

var_gamma <- sigma2* solve(t(x) %*% Z %*% W_opt %*% t(Z) %*% x)


stdvs_BGMM_sys  <-  sqrt(diag(var_gamma))
stdvs_BGMM_sys

tstats_BGMM_sys <- gamma/stdvs_BGMM_sys
tstats_BGMM_sys <- c(tstats_BGMM_sys)

n  <- length(y)
k  <- ncol(x)
df <- n-k

pvals_GGMM_sys   <- 2*(1-pt(abs(tstats_BGMM_sys),df))


#save output
names(gamma) <- colnames(x)
names(tstats_BGMM_sys) <- colnames(x)

#creating the table
coefs  <- round(gamma,3)
stdvs  <- round(stdvs_BGMM_sys,3)
tstats <- round(tstats_BGMM_sys,3)
pvals  <- round(pvals_GGMM_sys,3)


out_pvals_BGMM_sys = cbind(coefs,stdvs,tstats,pvals)

###################################################
############     Stacked Matrix     ###############
###################################################



N <-   length(unique(data$state))


#Let's loop over the states and create our big matrix

yit <- seq(1,T-2)
n_inst <- sum(yit)

#The chunks have size 28 !!!
#because yi27 is the last instrument of delta_ei29
y_1 <- data[1:(T-2),10]



#The number of columns of Z_1 is we would have without adding is 1 + 2 + .. + 27 = ninst
#but now for each chunk 1 + 2 + ... + 27, we had the 3 explanatory variables as instrument

Z_1 <- matrix(0,T-2,n_inst+((T-2)*3))
column <-2 

for (i in (1:(T-2))) {
  
  #Chunk is the bunch of yi that we are going to put in the Z_i matrix
  chunk <- as.numeric(y_1[1:i])
  #reverse chunk
  chunk <- rev(chunk)
  #print(chunk)
  exog_var <- data[i+2,16:18]

  
  endog_exog <- append(chunk,as.numeric(exog_var))
  print(endog_exog)
  
  for (j in (1:length(endog_exog))) {
    Z_1[i,(column+j-2)] <- endog_exog[j]
  }
  column <- column + i +3
}


x_1 <- exog_var

x_i <- x_1


Z_i <- matrix(0,T-2,n_inst+((T-2)*3))

Z <- Z_1


#We start at the column 30 and go by chunks of 29
#but we will only loop over the Y chunks until 27 !
#that is because yi27 will instrument delta_ei29
for (sst in seq(31, nrow(data), T)) {
  
  #print(sst)
  #The second chunk for example goes from 30 to 
  #30 + 27 
  y_i <- data[sst:(sst+27),10]
  
  x_i <- data[(sst+2):((sst+2)+28),16:18]
  
  #The matrix of instruments has size 27 x 378 ! 
  #each time period we use more lags
  Z_i <- matrix(0,T-2,n_inst+((T-2)*3))
  
  column <-2 
  #Next time I will start in row 57 which is row 59 -2 
  #and so one the keeps increase by 1 for each loop, this -2 takes the name iter
  
  
  for (i in (1:(T-2))) {
    
    #Chunk is the bunch of yi that we are going to put in the Z_i matrix
    chunk <- as.numeric(y_1[1:i])
    chunk <- rev(chunk)
    
    #print(chunk)
    exog_var <- x_i[i,]
    
    
    endog_exog <- append(chunk,as.numeric(exog_var))
    
    for (j in (1:length(endog_exog))) {
      Z_i[i,(column+j-2)] <- endog_exog[j]
    }
    column <- column + i +3
  }
  
  
  Z <- rbind(Z,Z_i)
  
}

#time to stack the matrix !

Z[Z == 0] <- NA

Z <- t(apply(Z,1,function(x){
  c(x[!is.na(x)],x[is.na(x)])}))


Z <- Z[,colSums(is.na(Z))<nrow(Z)]


Z[is.na(Z)] <- 0

n_inst_var <- N*(T-2)

#Imagine our big H is very big -> we try to create 
diagonal <- 2
offdiagonal<- -1
H <- matrix(0,n_inst_var,n_inst_var)
diag(H) <- diagonal
diag(H[-1,])<-offdiagonal
diag(H[,-1])<-offdiagonal

#this is a matrix with diagonal 2 and -1 on each side 

data_reg <- na.omit(data)

x1 <- as.matrix(na.omit(data_reg[,16:18]))
x2 <- as.matrix(na.omit(data_reg[,19]))
x<- cbind(x2,x1) #V1 like vendogemous

data_reg <- na.omit(data)
y <- as.matrix(na.omit(data_reg[,15]))

W_notinv <- t(Z) %*% H %*% Z 

#we have the big W optimal
W_opt <- solve(W_notinv,tol = 1e-22)


gamma <- solve(t(x) %*% Z %*% W_opt %*% t(Z) %*% x) %*%  t(x) %*% Z %*% W_opt %*% t(Z) %*% y

y_pred <- y
yhat   <- as.vector(x%*%gamma)
res    <- y_pred-yhat

sigma2 <- as.vector(t(res)%*%res)/df

sigma2 <- sigma2/2


var_gamma <- sigma2* solve(t(x) %*% Z %*% W_opt %*% t(Z) %*% x)




stdvs_BGMM_sys  <-  sqrt(diag(var_gamma))
stdvs_BGMM_sys

tstats_BGMM_sys <- gamma/stdvs_BGMM_sys
tstats_BGMM_sys <- c(tstats_BGMM_sys)

n  <- length(y)
k  <- ncol(x)
df <- n-k

pvals_GGMM_sys   <- 2*(1-pt(abs(tstats_BGMM_sys),df))


#save output
names(gamma) <- colnames(x)
names(tstats_BGMM_sys) <- colnames(x)

#creating the table
coefs  <- round(gamma,3)
stdvs  <- round(stdvs_BGMM_sys,3)
tstats <- round(tstats_BGMM_sys,3)
pvals  <- round(pvals_GGMM_sys,3)


out_pvals_BGMM_sys = cbind(coefs,stdvs,tstats,pvals)




