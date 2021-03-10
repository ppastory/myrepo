#Case 3 - GMM function

GMM_own = function(y,x,z,w)
  
# y is a (Nx1) vector of the dependent variable
# X is a (Nxk) vector of k explanatory variables 
# (for a constant include vector of 1's in the first column of NxK matrix
# NOTE: maybe add option in the function (line 3) to add a constant in the NxK matrix)
# Z is a (Nxr) vector of r instruments (r can be different than k)
# W is an option to select the optimal weighting matrix (0 for standard, 1 for White, 2 for Newey-West)
  
# The 1 step GMM estimator
  
{
  
  n <- length(y)
  k <- ncol(x)
  r <- ncol(z)
  df <- n-k 
  dr <- k-r
  
#define some variables
  
  kmoments = col_count(x)
  beta_vector = 
  if (beta_vector > kmoments) {
    print('GMM not identified')
  }
  
  
#count the moment conditions. Make sure that you have enough moment conditions
# to compute beta vector  
 
if (z==o && w==0) {
  
  betavector = inv(t(x)%*%x)%*%t(x)%*%y
  
  res = y - x%*%betavector
  var = t(res)%*%res%*%(n-K)

  
#compute the beta vector
  
  
}
  
  
   
xlarge = cbind(x,)
  coefs = inv(t(z)%*%x)%*%t(z)%*%y
  
  
  
    
}