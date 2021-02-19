####################################################################
##################               OLS              ##################
####################################################################
OLS_own = function (y,x) 
{
  n  <- length(y)
  k  <- ncol(x)
  df <- n-k
  
  ## Run OLS
  xy     <- t(x)%*%y
  xxi    <- solve(t(x)%*%x)
  coefs  <- as.vector(xxi%*%xy)

  yhat   <- as.vector(x%*%coefs)
  res    <- y-yhat
  sigma2 <- as.vector(t(res)%*%res/df)
  
  stdvs  <- sqrt(sigma2)*sqrt(diag(xxi))
  tstats <- coefs/stdvs
  pvals  <- 2*(1-pt(abs(tstats),df))

  ## Save output
  names(coefs) <- colnames(x)
  
  coefs  <- round(coefs,3)
  stdvs  <- round(stdvs,3)
  tstats <- round(tstats,3)
  pvals  <- round(pvals,3)
  
  out = rbind(coefs,stdvs,tstats,pvals)
  out = t(out)
  return(out)
}

