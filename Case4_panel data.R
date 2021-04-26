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
#Case 4.1. Monte Carlo simulation for OLS, X deterministic ########
####################################################################
# Read from file
data  <- read_excel("Data_Baltagi.xlsx")
Pdata <- pdata.frame(data, index=c("state", "year"))

#Estimate the dynamic model in levels using the FE estimator
#create lag of lnC_it
lag_lnC_it<-lag(Pdata$ln.C_it, k=1)
#remove year of 1963 because NA 



  