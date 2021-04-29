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
#install.packages("expm")


## Load required packages
library(plm) 
#library(expm)

## Load your own functions
source("Case1_Functions.R")
source("Case3_Functions.R")


####################################################################
###             Load data into the workspace                    ####
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

#ok now we've got all the data we need!

##########################################################################

#Compute the first difference model

FD_OLS = OLS_own(y_D_diff,x_D_diff,0)


