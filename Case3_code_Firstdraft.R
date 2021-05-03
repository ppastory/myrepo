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
if(!require(plm)){install.packages("plm")}
if(!require(dplyr)){install.packages("dplyr")}## panel data package
#install.packages("expm")


## Load required packages
library(plm)
library(dplyr)
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

y_D_diff = na.omit(y_D_diff) #Remove NA's otherwise the OLS functions returns an error
x_D_diff = na.omit(x_D_diff)

#ok now we've got all the data we need!

##########################################################################

#Compute the first difference model

FD_OLS = OLS_own(y_D_diff,x_D_diff,1)

#########################################################################

#Show that taking first differences automatically induces an endogeneity problem in this (first difference mode.red) model
#Answer see slide 73 of panel lecture. DeltaLNCit-1 and error term are mechanically correlated in a pooled OLS setting as is the model in line 87.

#########################################################################

#Which lagged values of the variables are valid instruments?

#construct lagged values
pdt$ln.C_itL2 <- lag(pdt$ln.C_it,2)
pdt$ln.C_itL3 <- lag(pdt$ln.C_it,3)
pdt$ln.C_itL4 <- lag(pdt$ln.C_it,4)
pdt$ln.C_itL5 <- lag(pdt$ln.C_it,5)
pdt$ln.P_itL1 <- lag(pdt$ln.P_it,1)
pdt$ln.P_itL2 <- lag(pdt$ln.P_it,2)
pdt$ln.P_itL3 <- lag(pdt$ln.P_it,3)
pdt$ln.P_itL4 <- lag(pdt$ln.P_it,4)
pdt$ln.P_itL5 <- lag(pdt$ln.P_it,5)
pdt$ln.Pn_itL1 <- lag(pdt$ln.Pn_it,1)
pdt$ln.Pn_itL2 <- lag(pdt$ln.Pn_it,2)
pdt$ln.Pn_itL3 <- lag(pdt$ln.Pn_it,3)
pdt$ln.Pn_itL4 <- lag(pdt$ln.Pn_it,4)
pdt$ln.Pn_itL5 <- lag(pdt$ln.Pn_it,5)
pdt$ln.Y_itL1 <- lag(pdt$ln.Y_it,1)
pdt$ln.Y_itL2 <- lag(pdt$ln.Y_it,2)
pdt$ln.Y_itL3 <- lag(pdt$ln.Y_it,3)
pdt$ln.Y_itL4 <- lag(pdt$ln.Y_it,4)
pdt$ln.Y_itL5 <- lag(pdt$ln.Y_it,5)

pdt_noNaN_instruments    <- na.omit(pdt)

#extract the dependent, exogenous, endogenous and instrumental variables
lnC_it_D_gmm    <- pdt_noNaN_instruments$ln.C_it
lnC_itL1_D_gmm  <- pdt_noNaN_instruments$ln.C_itL1
lnP_it_D_gmm    <- pdt_noNaN_instruments$ln.P_it
lnPn_it_D_gmm   <- pdt_noNaN_instruments$ln.Pn_it
lnY_it_D_gmm    <- pdt_noNaN_instruments$ln.Y_it

ln.C_itL2 <- pdt_noNaN_instruments$ln.C_itL2
ln.C_itL3 <- pdt_noNaN_instruments$ln.C_itL3
ln.C_itL4 <- pdt_noNaN_instruments$ln.C_itL4
ln.C_itL5 <- pdt_noNaN_instruments$ln.C_itL5
ln.P_itL1 <- pdt_noNaN_instruments$ln.P_itL1
ln.P_itL2 <- pdt_noNaN_instruments$ln.P_itL2
ln.P_itL3 <- pdt_noNaN_instruments$ln.P_itL3
ln.P_itL4 <- pdt_noNaN_instruments$ln.P_itL4
ln.P_itL5 <- pdt_noNaN_instruments$ln.P_itL5
ln.Pn_itL1 <- pdt_noNaN_instruments$ln.Pn_itL1
ln.Pn_itL2 <- pdt_noNaN_instruments$ln.Pn_itL2
ln.Pn_itL3 <- pdt_noNaN_instruments$ln.Pn_itL3
ln.Pn_itL4 <- pdt_noNaN_instruments$ln.Pn_itL4
ln.Pn_itL5 <- pdt_noNaN_instruments$ln.Pn_itL5
ln.Y_itL1 <- pdt_noNaN_instruments$ln.Y_itL1
ln.Y_itL2 <- pdt_noNaN_instruments$ln.Y_itL2
ln.Y_itL3 <- pdt_noNaN_instruments$ln.Y_itL3
ln.Y_itL4 <- pdt_noNaN_instruments$ln.Y_itL4
ln.Y_itL5 <- pdt_noNaN_instruments$ln.Y_itL5


#Group the extracted variables in vectors and matrices
y_D_gmm<-as.vector(lnC_it_D_gmm)
x_D_gmm<-as.matrix(cbind(Cnst=1,lnC_itL1_D_gmm,lnP_it_D_gmm,lnPn_it_D_gmm,lnY_it_D_gmm))
z_D_gmm<-as.matrix(cbind(ln.C_itL2,ln.C_itL3,ln.C_itL4,ln.C_itL5,ln.P_itL1,ln.P_itL2,ln.P_itL3,ln.P_itL4,ln.P_itL5,ln.Pn_itL1,ln.Pn_itL2,ln.Pn_itL3,ln.Pn_itL4,ln.Pn_itL5,ln.Y_itL1,ln.Y_itL2,ln.Y_itL3,ln.Y_itL4,ln.Y_itL5))

# Take first differences, no need to difference the lags (see case slide 10)
lnC_it_D_diff_gmm   = diff(lnC_it_D, differences = 1)
lnC_itL1_D_diff_gmm = diff(lnC_itL1_D, differences = 1)
lnP_it_D_diff_gmm   = diff(lnP_it_D, differences = 1)                                    # Take first differences                
lnPn_it_D_diff_gmm  = diff(lnPn_it_D, differences = 1)
lnY_it_D_diff_gmm   = diff(lnY_it_D, differences = 1)

y_D_diff_gmm<-as.vector(lnC_it_D_diff_gmm)
x_D_diff_gmm<-as.matrix(cbind(Cnst=1,lnC_itL1_D_diff_gmm,lnP_it_D_diff_gmm,lnPn_it_D_diff_gmm,lnY_it_D_diff_gmm))

y_D_diff_gmm = na.omit(y_D_diff_gmm) #Remove NA's otherwise the OLS functions returns an error
x_D_diff_gmm = na.omit(x_D_diff_gmm)

#check correlation of instruments with lnC_itL1_D_diff_gmm (the endogenous variable)
lnC_it_D_diff_gmm_NoNaN = na.omit(lnC_itL1_D_diff_gmm)
test <- subset.data.frame(pdt_noNaN, year==67 | year==68|year==69|year==70|year==71|year==72|year==73|year==74|year==75|year==76|year==77|year==78|year==79|year==80|year==81|year==82|year==83|year==84|year==85|year==86|year==87|year==88|year==89|year==90|year==91|year==92)
lnC_itL1_D_test  <- test$ln.C_itL1
lnC_itL1_D_diff_test = diff(lnC_itL1_D_test, differences = 1)
lnC_itL1_D_diff_test = as.vector(lnC_itL1_D_diff_test)
lnC_itL1_D_diff_test = na.omit(lnC_itL1_D_diff_test)


Covxz <- cov(lnC_itL1_D_diff_test,z_D_gmm)
Corrxz <- cor(lnC_itL1_D_diff_test, z_D_gmm)

#check correlation with residuals FD_OLS
lnC_it_D_diff_gmm_NoNaN = na.omit(lnC_it_D_diff_gmm)
lnC_it_D_diff_gmm_NoNaN_subset <- subset.data.frame(lnC_it_D_diff_gmm_NoNaN, year == 68)


pdt2 <- pdata.frame(data)
pdt2 <- transform(pdt2, ln.C_itL1_subset=lag(lnC_it,1))
pdt2 <- transform(pdt2, ln.C_it_D_diff_subset=diff(ln.C_it,1))
pdt2 <- transform(pdt2, ln.P_it_D_diff_subset=diff(ln.P_it,1))
pdt2 <- transform(pdt2, ln.Pn_it_D_diff_subset=diff(ln.Pn_it,1))
pdt2 <- transform(pdt2, ln.Y_it_D_diff_subset=diff(ln.Y_it,1))

pdt2_subset <- subset.data.frame(pdt2, year==68|year==69|year==70|year==71|year==72|year==73|year==74|year==75|year==76|year==77|year==78|year==79|year==80|year==81|year==82|year==83|year==84|year==85|year==86|year==87|year==88|year==89|year==90|year==91|year==92)

#extract first differenced dependent and independent variables from the subset of pdt (pdt2_subset)

lnC_it_D_diff_subset    <- pdt2_subset$ln.C_it_D_diff_subset
lnC_itL1_D_diff_subset  <- pdt2_subset$ln.C_itL1_subset
lnP_it_D_diff_subset    <- pdt2_subset$ln.P_it_D_diff_subset
lnPn_it_D_diff_subset   <- pdt2_subset$ln.Pn_it_D_diff_subset
lnY_it_D_diff_subset   <- pdt2_subset$ln.Y_it_D_diff_subset

y_D_diff_subset<-as.vector(lnC_it_D_diff_subset)
x_D_diff_subset<-as.matrix(cbind(Cnst=1,lnC_itL1_D_diff_subset,lnP_it_D_diff_subset,lnPn_it_D_diff_subset,lnY_it_D_diff_subset))


FD_OLS_test <- OLS_own(y_D_diff_subset,x_D_diff_subset,0) # I cannot retrieve White residuals in case 1 code #also error in ols formula when regressing two vectors
FD_OLS_test$residuals
FD_OLS_test$estimation
residuals_FD_OLS <- FD_OLS_test$residuals

Corrzres = cor(residuals_FD_OLS, z_D_gmm)

#check correlation matrix with other exogenous variables
correxz = cor(x_D_diff_subset[,2:5],z_D_gmm)

#the second lag of consumption, third lag of price and third lag of Pn are candidates of lagged exogenous variables 
#1) correlated with endogenous variable, uncorrelated with error term, low correlation with the remaining exogenous variables

#Test validity lags with Stock and Yogo, apply to the exactly identified model
Z <- cbind(x_D_diff_subset[,3:5],z_D_gmm[,1])
StockYogo_secondlag_lnCit = OLS_own(lnC_itL1_D_diff_subset, Z,0)
Z <- cbind(x_D_diff_subset[,3:5],z_D_gmm[,7])
StockYogo_thirdlag_lnPit = OLS_own(lnC_itL1_D_diff_subset, cbind(x_D_diff_subset[,3:5],Z),0)
Z <- cbind(x_D_diff_subset[,3:5],z_D_gmm[,12])
StockYogo_thirdlag_lnPnit = OLS_own(lnC_itL1_D_diff_subset, cbind(x_D_diff_subset[,3:5],Z),0)

#compute mu2 of Stock and Yogo

alpha <- StockYogo_secondlag_lnCit$estimation[,1]
Z <- cbind(x_D_diff_subset[,3:5],z_D_gmm[,1])
R <- ncol(Z)
muepsilon2 <- var(StockYogo_secondlag_lnCit$residuals)

Fstatistic1 = t(alpha)%*%t(Z)%*%Z%*%alpha/R/muepsilon2

alpha <- StockYogo_thirdlag_lnPit$estimation[,1]
Z <- cbind(x_D_diff_subset[,3:5],z_D_gmm[,7])
R <- ncol(Z)
muepsilon2 <- var(StockYogo_thirdlag_lnPit$residuals)

Fstatistic2 = t(alpha)%*%t(Z)%*%Z%*%alpha/R/muepsilon2

alpha <- StockYogo_thirdlag_lnPnit$estimation[,1]
Z <- cbind(x_D_diff_subset[,3:5],z_D_gmm[,12])
R <- ncol(Z)
muepsilon2 <- var(StockYogo_thirdlag_lnPnit$residuals)

Fstatistic3 = t(alpha)%*%t(Z)%*%Z%*%alpha/R/muepsilon2

if(Fstatistic1>10)
{print("strong instrument")} else 
{print("weak instrument")}

#estimate exactly identified model using the second lag of lnCit as instrument
GMM_exact = GMM_own(y_D_diff_subset, x_D_diff_subset, z_D_gmm[,1],0)

#estimate overidentified model
z_over = cbind(z_D_gmm[,1], z_D_gmm[,7], z_D_gmm[,12])
GMM_over = GMM_own(y_D_diff_subset, x_D_diff_subset, z_D_gmm[,1],1)


#in exactly identified model it is not possible to test the validity of the moment conditions as the moment conditions are
#identifying. In the overidentified model it is possible to test for the validity of the moment conditions using the hansen test. 





