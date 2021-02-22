####################################################################
##################               OLS              ##################
####################################################################

OLS_own = function (y,x,W) 
#W is the option for different corrections
#1 for White
#2 for GLS

{
  n  <- length(y)
  k  <- ncol(x)
  df <- n-k
  W  <- 1:2
  
  ## Run OLS
  xy     <- t(x)%*%y
  xxi    <- solve(t(x)%*%x)
  coefs  <- as.vector(xxi%*%xy)

  yhat   <- as.vector(x%*%coefs)
  res    <- as.vector(y-yhat)
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
####################################################################
###############     OLS with White correction     ##################
####################################################################

if W=1
{
  res2       <- as.vector(res%*%res)
  sigma2omega<- diag(res2)
  cov_W      <- xxi * t(x) * sigma2omega * x * xxi
  stdvs_w    <- sqrt(diag(cov_W))
  tstats_w   <- coefs/stdvs_w
  pvals_w    <- 2*(1-pt(abs(tstats_w),df))
}


####################################################################
##################               GLS              ##################
####################################################################

if W=2
{
  P         <- diag(sqrt(res*res)^(-1))
  omega     <- P^(-1)*t(P)^(-1)
  coefs_GLS <- (t(x)*omega^(-1)*x)^(-1)*t(x)*omega^(-1)*y
}

