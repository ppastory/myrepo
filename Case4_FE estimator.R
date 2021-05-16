####################################################################
##################              FE estimator             ##################
####################################################################

FE = function (data,TFE) 
###build matrix D
{
  y <- as.matrix(data[,10])
  
  x <- as.matrix(data[,11:14])
  
  
  idx <- sort(unique(data$state))
  state_i <- matrix(0, nrow = nrow(data), ncol = length(idx))
  
  for (j in 1:length(idx)) { 
    state_i[,j] <- as.integer(data$state == idx[j])
  }
  
  idx <- sort(unique(data$year))
  year_i <- matrix(0, nrow = nrow(data), ncol = length(idx))
  
  for (j in 1:length(idx)) { 
    year_i[,j] <- as.integer(data$year == idx[j])
  }
  
  #Let's write the demeaning matrix by individuals
  #I need the identity matrix
  NT <- nrow(data)
  I   <- matrix(0, NT, NT)
  diag(I) <- 1
  
  P_i <- state_i %*% solve(t(state_i) %*% state_i) %*% t(state_i)
  
  P_T <- year_i %*% solve(t(year_i) %*% year_i) %*% t(year_i)
  
  y_demean <- y - P_i %*% y
  
  x_demean <- x - P_i %*% x
  
  
  y_demean_it <- y_demean - P_T %*% y_demean
  
  x_demain_it <- x_demean - P_T %*% x_demean
  
  if (TFE == 1){
    x <- x_demain_it
    y <- y_demean_it
  } else{
    x <- x_demean
    y <- y_demean
    
  }

  xy     <- t(x)%*%y
  xxi    <- solve(t(x)%*%x) #this is (X' X)^(-1)
  coefs  <- as.vector(xxi%*%xy)
  
  yhat   <- as.vector(x%*%coefs)
  res    <- y-yhat
  
  sigma2 <- as.vector(t(res)%*%res)/df
  
  #case where we only want to do 
  
  var <- sigma2*xxi
  stdvs <- sqrt(diag(var))
  
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
  
  output <- list("estimation" = out, "data" = cbind(y,x))
  
  return(output)
  
}