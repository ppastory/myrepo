####################################################################
##################               OLS              ##################
####################################################################

OLS_own = function (y,x,w) 
#W is the option for different corrections
#1 for White
#2 for GLS

{
  n  <- length(y)
  k  <- ncol(x)
  df <- n-k
  
  n  <- length(y)
  k  <- ncol(x)
  df <- n-k
  
  ## Run OLS
  xy     <- t(x)%*%y
  xxi    <- solve(t(x)%*%x) #this is (X' X)^(-1)
  coefs  <- as.vector(xxi%*%xy)
  
  yhat   <- as.vector(x%*%coefs)
  res    <- y-yhat
  sigma2 <- as.vector(t(res)%*%res/df)
  
  print (w)
  #case where we only want to do 
  
  if (w == 0){
    
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
    
  } else if (w == 1){  #this is the case with white correction
    
    res2       <- res%*%t(res)
    #the non-diagnonal elemets of res2 need to have 0
    #I take away the diagonals
    diagonal <- diag(res2) 
    #I create a matrix of 0 of diam
    res2 <- matrix(0,nrow(res2), ncol(res2)) 
    #I put back the diagonal element in the diagnonals
    diag(res2) <- diagonal
    #res2 is sigma2 omega in the formulas
    cov_W      <- xxi %*% t(x) %*% res2 %*% x %*% xxi
    cov_W
    stdvs_w    <- sqrt(diag(cov_W))
    stdvs_w
    tstats_w   <- coefs/stdvs_w
    tstats_w
    pvals_w    <- 2*(1-pt(abs(tstats_w),df))
    pvals_w
    #save output
    names(coefs) <- colnames(x)
    
    coefs  <- round(coefs,3)
    stdvs  <- round(stdvs_w,3)
    tstats <- round(tstats_w,3) #this vector is a row and it should be a column
    pvals  <- round(pvals_w,3)
    
    out = rbind(coefs,stdvs,tstats,pvals)
    out = t(out)
    out
    return(out)
    
  } else if (w == 2){  #this is the case with GLS
    
    res2       <- res%*%t(res)
    #the non-diagnonal elemets of res2 need to have 0
    #I take away the diagonals
    diagonal <- diag(res2) 
    #I create a matrix of 0 of diam
    res2 <- matrix(0,nrow(res2), ncol(res2)) 
    #I put back the diagonal element in the diagnonals
    #but this time I take the inverse to get P
    diag(res2) <- 1/sqrt(diagonal)
    #res2 is sigma2 omega in the formulas
    P <- res2
    
    omega_1    <- t(P)%*%P
    
    
    coefs_GLS <- solve((t(x) %*% omega_1 %*% x)) %*% t(x) %*% omega_1 %*% y
    coefs_GLS <- c(coefs_GLS)
    
    cov_GLS   <- sigma2 * solve(t(x) %*% omega_1 %*%x)
    cov_GLS
    
    stdvs_GLS    <- sqrt(diag(cov_GLS))
    stdvs_GLS
    
    tstats_GLS   <- coefs_GLS/stdvs_GLS
    tstats_GLS <- c(tstats_GLS)
    
    pvals_GLS    <- 2*(1-pt(abs(tstats_GLS),df))
    
    
    #save output
    names(coefs_GLS) <- colnames(x)
    names(tstats_GLS) <- colnames(x)
    
    #creating the table
    coefs  <- round(coefs_GLS,3)
    stdvs  <- round(stdvs_GLS,3)
    tstats <- round(tstats_GLS,3)
    pvals  <- round(pvals_GLS,3)
    
    
    out_GLS = rbind(coefs,stdvs,tstats,pvals)
    
    out_GLS = t(out_GLS)
    return(out_GLS)
    
  } else if (w == 3){  #this is the case with EGLS
    
    res2       <- res%*%t(res)
    #I want a proxy for the errors squares 
    diagonal <- diag(res2) 
    
    ln_sigma <- log(diagonal)
    ln_x <- x #Do we have to log the x if it is already done?
    
    xy     <- t(ln_x)%*%ln_sigma
    xxi    <- solve(t(ln_x)%*%ln_x) #this is (X' X)^(-1)
    coefs  <- as.vector(xxi%*%xy)
    
    yhat   <- as.vector(ln_x%*%coefs)
    sigma_hat <- exp(yhat)
    
    #I want to put my estimated sigma back in a white correction matrix
    res2 <- matrix(0,nrow(res2), ncol(res2)) 
    #I put back the diagonal element in the diagnonals
    #but this time I take the inverse to get P
    diag(res2) <- sigma_hat
    
    cov_EGLS      <- xxi %*% t(x) %*% res2 %*% x %*% xxi
    cov_EGLS
    
    coefs_EGLS <- solve((t(x) %*% res2 %*% x)) %*% t(x) %*% res2 %*% y
    coefs_EGLS <- c(coefs_GLS)
    
    stdvs_EGLS   <- sqrt(diag(cov_W))
    stdvs_EGLS
    
    tstats_EGLS   <- coefs/stdvs_w
    tstats_EGLS
    
    pvals_EGLS    <- 2*(1-pt(abs(tstats_w),df))
    pvals_EGLS
    
    #save output
    #names(coefs) <- colnames(x)
    
    #creating the table
    coefs  <- round(coefs_EGLS,3)
    stdvs  <- round(stdvs_EGLS,3)
    tstats <- round(tstats_EGLS,3)
    pvals  <- round(pvals_EGLS,3)
    
    
    out_EGLS = rbind(coefs,stdvs,tstats,pvals)
    
    out_EGLS = t(out_EGLS)
    return(out_EGLS)
    
  }
  
  
}



