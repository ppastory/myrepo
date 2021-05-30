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
#Case 4.2. GMMd                                       ########
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