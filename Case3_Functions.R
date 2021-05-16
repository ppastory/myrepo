#Case 3 - GMM function

GMM_own = function(y,x,z,w)
  
# y is a (Nx1) vector of the dependent variable
# X is a (Nxk) vector of k explanatory variables 
# (for a constant include vector of 1's in the first column of NxK matrix
# NOTE: maybe add option in the function (line 3) to add a constant in the NxK matrix)
# Z is a (Nxr) vector of r instruments (r can be different than k). Including exogenous explanatory variables.
# W is an option to select the optimal weighting matrix (0 for standard, 1 for White, 2 for Newey-West)
  
# The 1 step GMM estimator
  
{
  
  n <- length(y)
  k <- ncol(x)
  r <- ncol(z)
  df <- n-k 
  dr <- k-r
  
#count the moment conditions. Make sure that you have enough moment conditions
#to compute betahat
 
#weighting matrix only for overidentified
  if (w==-1) {
  
    betahat = solve(t(x)%*%x)%*%t(x)%*%y
    
    reshat = y - x%*%betahat
    
    sigmahat2 = t(reshat)%*%reshat%*%(n-K)
    
    xxi = solve(t(x)%*%x%*%solve(t(x)%*%x)%*%t(x)%*%x)
    
    covarbeta = sigmahat2 * xxi
  
    stderror < - sqrt(diag(covarbeta))
    
    tstats <- betahat/stderror
    
    pvals <- 2*(1-pt(abs(tstats),df)) 
    
    #Note test this again with t.test function in R. I do not understand the notation.
    
    ## Save output
    names(betahat) <- colnames(x)
    
    betahat  <- round(betahat,3)
    stderror  <- round(stderror,3)
    tstats <- round(tstats,3)
    pvals  <- round(pvals,3)
    
    out = rbind(betahat,stderror,tstats,pvals)
    out = t(out)
    
    output <- list("estimation" = out, "residuals" = res, "param" = k)
    
    return(output)
  
    }
    else if (w==0) {
    
    #Note: build a check such that x and z have the same dimensions
    #assume exactly identified model so X and K matrix have the same dimensions
    #GIV also takes care of exactly identified see sl. 57
 
    wmatrix = solve(t(z)%*%z) #when errors are iid 
    
    betahativ = solve(t(x)%*%z%*%wmatrix%*%t(z)%*%x)%*%t(x)%*%z%*%solve(t(z)%*%z)%*%t(z)%*%y
    
    betahativ <- as.vector(betahativ) 
    
    #compute GMM finite sample standard errors
    res = y - x%*%betahativ
    
    sigmahat2 = as.vector((t(res)%*%res))/(n-k)
    
    xxi = solve(t(x)%*%z%*%solve(t(z)%*%z)%*%t(z)%*%x)
    
    covarbetaiv = sigmahat2 * xxi
    
    stderror <- sqrt(diag(covarbetaiv))
    
    tstats <- betahativ/stderror
    
    pvals <- 2*(1-pt(abs(tstats),df)) 
    
    ## Save output
    names(betahativ) <- colnames(x)
    
    betahativ  <- round(betahativ,3)
    stderror  <- round(stderror,3)
    tstats <- round(tstats,3)
    pvals  <- round(pvals,3)
    
    out = rbind(betahativ,stderror,tstats,pvals)
    out = t(out)
    
    output <- list("estimation" = out, "residuals" = res, "param" = k)
    
  return(output)
    }  
  else if (w==1) {
    
  #for one step GMM is always for iid errors, when not errors are not iid then choose for two step. 
  
  wmatrix = solve(t(z)%*%z) #to compute sigma2omega
  
  betahativ = solve(t(x)%*%z%*%wmatrix%*%t(z)%*%x)%*%t(x)%*%z%*%solve(t(z)%*%z)%*%t(z)%*%y
  
  betahativ <- as.vector(betahativ) 
  
  res = y - x%*%betahativ
  
  res2 = res%*%t(res)
  
  diagonal <- diag(res2)
  
  res2 = matrix(0,nrow(res2), ncol(res2))
  
  diag(res2) <- diagonal
  
  sigma2omega = res2
          
  betahativ = solve(t(x)%*%z%*%solve(t(z)%*%sigma2omega%*%z)%*%t(z)%*%x)%*%t(x)%*%z%*%solve(t(z)%*%sigma2omega%*%z)%*%t(z)%*%y
  
  betahativ <- as.vector(betahativ) 
  
  #compute GMM standard errors with White weighting matrix
  
  res = y - x%*%betahativ
  
  sigmahat2 = as.vector((t(res)%*%res))/(n-k) #check if k is correct? Because we can have more than k parameters
  
  xxi = solve(t(x)%*%z%*%solve(t(z)%*%sigma2omega%*%z)%*%t(z)%*%x)
  
  covarbetaiv = sigmahat2*xxi
  
  #Compute Windmeijer correction (nice to have)
  
  stderror <- sqrt(diag(covarbetaiv))
  
  tstats <- betahativ/stderror
  
  pvals <- 2*(1-pt(abs(tstats),df)) 
  
  #Compute Hansen diagnostic test for identified models
  
  if (r>k) {
  
  J1 = t(t(z)%*%res)%*%solve(t(z)%*%sigma2omega%*%z)%*%(t(z)%*%res)
  
  #define z2 for a subset of instruments Check this line
  
  J2 = t(t(z2)%*%res)%*%solve(t(z2)%*%sigma2omega%*%z2)%*%(t(z2)%*%res)
  
  diffJ= J1 - J2
  
  #compute number of suspect moment conditions
  
  degreesoffreedomchisquared = r-f
  
  #test joint validity of moment conditions using J stats
  
  pval = dchisq(diffJ,degreesoffreedomchisquared)
  
  }
  
  ## Save output
  names(betahativ) <- colnames(x)
  
  betahativ  <- round(betahativ,3)
  stderror  <- round(stderror,3)
  tstats <- round(tstats,3)
  pvals  <- round(pvals,3)
  
  out = rbind(betahativ,stderror,tstats,pvals)
  out = t(out)
  
  output <- list("estimation" = out, "residuals" = res, "param" = k)
  
  return(output)  
  }
  else if (w==2) {
  
    
    #for one step GMM is always for iid errors, when not errors are not iid then choose for two step. 
    
    #compute Newey-West Weighting matrix i.e. HAC standard errors
    
    wmatrix = solve(t(z)%*%z) #to compute sigma2omega
    
    betahativ = solve(t(x)%*%z%*%wmatrix%*%t(z)%*%x)%*%t(x)%*%z%*%solve(t(z)%*%z)%*%t(z)%*%y
    
    betahativ <- as.vector(betahativ) 
    
    res = y - x%*%betahativ
    
    res2 = res%*%t(res)
    
    diagonal <- diag(res2)
    
    res2 = matrix(0,nrow(res2), ncol(res2))
    
    diag(res2) <- diagonal
    
    #sigma2omega = res2
    
    #Define Newey-West lag
    L = n^(1/4)
    
    #Define wl for Newey-West
    for (i in 1:n){
      muxx = res2[i,i]%*%z[i,]%*%t(z[,i])
    }
    
    for (l in 1:L){
      for (j in (l+1):n){
        wl = 1 - (l/(L+1))
        doublesum = wl%*%res[j]%*%res[l]%*%(z[t]%*%t(z[,j-l])+z[j-l]%*%t(z[,j]))
      }
    }
    
    neweywestcovar = muxx + doublsum
    
    betahativ = solve(t(x)%*%z%*%solve(t(z)%*%neweywestcovar%*%z)%*%t(z)%*%x)%*%t(x)%*%z%*%solve(t(z)%*%neweywestcovar%*%z)%*%t(z)%*%y
    betahativ <- as.vector(betahativ) 
    
    #compute GMM standard errors with White weighting matrix
    
    res = y - x%*%betahativ
    
    sigmahat2 = as.vector((t(res)%*%res))/(n-k) #check if k is correct? Because we can have more than k parameters
    
    xxi = solve(t(x)%*%z%*%solve(t(z)%*%neweywestcovar%*%z)%*%t(z)%*%x)
    
    covarbetaiv = sigmahat2%*%xxi
    
    #Compute Windmeijer correction (nice to have)
    
    stderror <- sqrt(diag(covarbetaiv))
    
    tstats <- betahativ/stderror
    
    pvals <- 2*(1-pt(abs(tstats),df))
    
    #Compute Hansen diagnostic test
    
    if (r>k) {
      
      J1 = t(t(z)%*%res)%*%solve(t(z)%*%sigma2omega%*%z)%*%(t(z)%*%res)
      
      #define z2 for a subset of instruments
      
      J2 = t(t(z2)%*%res)%*%solve(t(z2)%*%sigma2omega%*%z2)%*%(t(z2)%*%res)
      
      diffJ= J1 - J2
      
      #compute number of suspect moment conditions
      
      degreesoffreedomchisquared = r-f
      
      #test joint validity of moment conditions using J stats
      
      pval = dchisq(diffJ,degreesoffreedomchisquared)
      
    }
    
    ## Save output
    names(betahativ) <- colnames(x)
    
    betahativ  <- round(betahativ,3)
    stderror  <- round(stderror,3)
    tstats <- round(tstats,3)
    pvals  <- round(pvals,3)
    
    out = rbind(betahativ,stderror,tstats,pvals)
    out = t(out)
    
    output <- list("estimation" = out, "residuals" = res, "param" = k)
    
    return(output) 
    
    
  }
    
}
