## Sample R script 
## Using the pooled OLS estimator to estimate the Cigarette demand example from the Baltagi textbook

#version 1.2

####################################################################
###              Set-up for use                                  ###
####################################################################
rm(list = ls())   # Clear workspace

## Set input file directory
setwd("/Users/geeverae/dropbox/My_doc/Onderwijs/Advanced Econometrics/R/Sample R code") #Mac
#setwd("C:\\Users\\...") #Windows

## Set output file directory: change to your own directory
output="/Users/geeverae/dropbox/My_doc/Onderwijs/Advanced Econometrics/R/Sample R code" #Mac
#output = "C:\\Users\\..." #Windows

## Install required packages (code is only executed if the packages are not yet installed)
if(!require(plm)){install.packages("plm")} ## panel data package

## Load required packages
library(plm)  

## Load your own functions
source("OwnFunctions.R")


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

####################################################################
###                Run OLS on static model                       ###
####################################################################

OLS_static = OLS_own(y_S,x_S)


####################################################################
###                Run OLS on static model                       ###
####################################################################

OLS_dynamic = OLS_own(y_D,x_D)
