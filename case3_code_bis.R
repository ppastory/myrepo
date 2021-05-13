####################################################################
###              Case 3 code                                     ###
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
setwd("C:\\Users\\qiszhang\\Documents\\programming/myrepo") #Windows

## Set output file directory: change to your own directory
#I want the output directory to be the same 
output = getwd() 
#so now everything is in one directory

## Install required packages (code is only executed if the packages are not yet installed)
if(!require(plm)){install.packages("plm")} ## panel data package
install.packages("data.table")


## Load required packages
library(plm) 
library(data.table)

## Load your own functions
source("Case1_Functions.R")
source("Case3_Functions.R")


####################################################################
###             Load data into the workspace                    ####
####################################################################
# Read from file
data    = read.table("Data_Baltagi.csv", header = TRUE,sep = ";")

colnames(data)[1] <- c('state')


#now let's create the lags :) 

lg1 <- function(x)c(NA, x[1:(length(x)-1)])

data$ln.C_it_1 <- unlist(tapply(data$`ln.C_it`, data$state, lg1))
data<-na.omit(data)
data$ln.C_it_2 <- unlist(tapply(data$`ln.C_it_1`, data$state, lg1))
data<-na.omit(data)
data$ln.C_it_3 <- unlist(tapply(data$`ln.C_it_2`, data$state, lg1))
data<-na.omit(data)
data$ln.C_it_4 <- unlist(tapply(data$`ln.C_it_3`, data$state, lg1))
data<-na.omit(data)
data$ln.C_it_5 <- unlist(tapply(data$`ln.C_it_4`, data$state, lg1))
data<-na.omit(data)

##First let's first difference our data!

data <- transform(data, dlnC_it = ave(`ln.C_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnP_it = ave(`ln.P_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnPn_it = ave(`ln.Pn_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnY_it = ave(`ln.Y_it`, state, FUN = function(x) c(NA, diff(x))))
data <-   transform(data, dlnC_it_1 = ave(ln.C_it_1, state, FUN = function(x) c(NA, diff(x))))

data<-na.omit(data)



######################
# validity condition #
######################

#we need to check if they are valid (exogenous)

y <- data[,12]

x <- as.matrix(data[,c(13:16)])

OLS_fd <- OLS_own(y,x,0)

res <- OLS_fd$residuals

cor(res,data[,c(7:11)])

#the correlation is not super close to 0 but let's have a more rigorous approach


######################
# relevance condition#
######################


#Let's remove the column we don't need
data <- data[-c(3:9)]

#let's see the correlation of our endogenous variable with the lags 
cor(data[,ncol(data)],data[,c(7:11)])

#not 0 they are relevant

#the correlation are small, let's have a more rigorous approach

##we need two instruments

#the exactly identified
#we use as instrument the first lag since it is the one correlating the most

Z_ei <- cbind(data[,c(13:16,7)])

Z_oi <- as.matrix(cbind(data[,c(13:16,7:11)]))


#we will do a F-test on our set of instruments (sl32)

x_endog <- data[,ncol(data)]


OLS_F1 <- OLS_own(x_endog, Z_oi,0)

OLS_F1



