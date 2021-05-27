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
source("Case4_FE Function.R")
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

FE_est$estimation

