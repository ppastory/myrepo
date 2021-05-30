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
 

  if (w==0) {
    
    
    #Note: build a check such that x and z have the same dimensions
    #assume exactly identified model so X and K matrix have the same dimensions
    #GIV also takes care of exactly identified see sl. 57
 
    wmatrix = solve(t(z)%*%z) #when errors are iid 
    
    betahativ = solve(t(x)%*%z%*%wmatrix%*%t(z)%*%x)%*%t(x)%*%z%*%wmatrix%*%t(z)%*%y
    
    betahativ <- as.vector(betahativ) 
    
    #compute GMM finite sample standard errors
    res = y - x%*%betahativ
    
    sigmahat2 = as.vector((t(res)%*%res))/(n-k)
    
    xxi = solve(t(x)%*%z%*%wmatrix%*%t(z)%*%x)
    
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
    
  #One step GMM is always for iid errors, when errors are not iid then choose for two step. 
  
  wmatrix = solve(t(z)%*%z) #to compute sigma2omega
  
  betahativ = solve(t(x)%*%z%*%wmatrix%*%t(z)%*%x)%*%t(x)%*%z%*%wmatrix%*%t(z)%*%y
  
  betahativ <- as.vector(betahativ) 
  
  res = y - x%*%betahativ
  
  res2 = res%*%t(res)
  
  diagonal <- diag(res2)
  
  res2 = matrix(0,nrow(res2), ncol(res2))
  
  diag(res2) <- diagonal
  
  #Define sigma2omega as component of Wmatrixhat
  sigma2omega = res2 
          
  betahativ = solve(t(x)%*%z%*%solve(t(z)%*%sigma2omega%*%z)%*%t(z)%*%x)%*%t(x)%*%z%*%solve(t(z)%*%sigma2omega%*%z)%*%t(z)%*%y
  
  betahativ <- as.vector(betahativ) 
  
  #compute GMM standard errors with 'White' weighting matrix
  
  res = y - as.vector(x%*%betahativ)
  
  sigmahat2 = as.vector((t(res)%*%res))/(n-k)
  
  xxi = solve(t(x)%*%z%*%solve(t(z)%*%sigma2omega%*%z)%*%t(z)%*%x)
  
  covarbetaiv = sigmahat2*xxi
  
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
  
  #Compute Sargan diagnostic test
  m <- t(z)%*%z
  
  w <- var(res) * m
  
  J = t(t(z)%*%res)%*%solve(w)%*%(t(z)%*%res)
  
  dfchisquared = r-k
  
  pvalj = dchisq(J,dfchisquared)
  
  output <- list("estimation" = out, "residuals" = res, "sargan test" = pvalj)
  
  return(output)
  
  }
  else if (w==2) {
  
    
    #for one step GMM is always for iid errors, when errors are not iid then choose for two step. 
    
    #compute Newey-West Weighting matrix i.e. HAC standard errors
    
    wmatrix = solve(t(z)%*%z) #To obtain initial estimate of sigma2
    
    betahativ = solve(t(x)%*%z%*%wmatrix%*%t(z)%*%x)%*%t(x)%*%z%*%wmatrix%*%t(z)%*%y
    
    betahativ <- as.vector(betahativ) 
    
    res = y - x%*%betahativ
    
    res2 = res*res

        
    n=length(z[,1])
    
    #Set up the Newey-West estimator. Define Newey-West lag
    L = n^(1/4)
    
    muxx = matrix(0,ncol(z),ncol(z))
    
    doublesum = matrix(0,ncol(z),ncol(z)) 
    
    #Define wl parameter for Newey-West
    for (i in 1:n){
      muxx_temporary = res2[i]*(z[i,]%*%t(z[i,]))
      muxx = muxx + muxx_temporary
    }
    
    for (l in 1:L){
      for (j in (l+1):n){
        wl = 1 - (l/(L+1))
        doublesum_temporary = wl*res[j]*res[j-l]*((z[j,]%*%t(z[j-l,]) + z[j-l,]%*%t(z[j,])))
        doublesum = doublesum + doublesum_temporary
      }
    }
    
    #Compute Newey-West estimator
    neweywestcovar = muxx + doublesum
    
    #Compute the second step GMM beta slide 67
    betahativ = solve(t(x)%*%z%*%solve(neweywestcovar)%*%t(z)%*%x)%*%t(x)%*%z%*%solve(neweywestcovar)%*%t(z)%*%y
    betahativ <- as.vector(betahativ)
    
    #compute GMM standard errors with Newey-West weighting matrix
    
    res = y - as.vector(x%*%betahativ)
    
    sigmahat2 = as.vector((t(res)%*%res))/(n-k)
    
    xxi = solve(t(x)%*%z%*%solve(neweywestcovar)%*%t(z)%*%x)
    
    covarbetaiv = sigmahat2*xxi
    
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
    
    
    
    #Compute Sargan diagnostic test
    m <- t(z)%*%z
    
    w <- var(res) * m
    
    J = t(t(z)%*%res)%*%solve(w)%*%(t(z)%*%res)
    
    dfchisquared = r-k
    
    pvalj = dchisq(J,dfchisquared)
    
    output <- list("estimation" = out, "residuals" = res, "sargan test" = pvalj)
    return(output)
    
    }
    
    
  }
    

