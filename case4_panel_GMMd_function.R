dGMM = function (data,full_set,stack,max_lag){
  #full set = 1 -> you 
  
  
  
  if (full_set == 1){
    #Extract T and N
    
    N <-   length(unique(data$state))
    T <-   length(unique(data$year)) #the real T
    
    
    #Let's loop over the states and create our big matrix
    
    yit <- seq(1,T-2)
    n_inst <- sum(yit)
    
    #The chunks have size 28 !!!
    #because yi27 is the last instrument of delta_ei29
    y_1 <- data[1:(T-2),10]
    
    
    
    #The number of columns of Z_1 is we would have without adding is 1 + 2 + .. + 27 = ninst
    #but now for each chunk 1 + 2 + ... + 27, we had the 3 explanatory variables as instrument
    
    Z_1 <- matrix(0,T-2,n_inst+((T-2)*3))
    column <-2 
    
    for (i in (1:(T-2))) {
      
      #Chunk is the bunch of yi that we are going to put in the Z_i matrix
      chunk <- as.numeric(y_1[1:i])
      #print(chunk)
      exog_var <- data[i+2,16:18]
      
      endog_exog <- append(as.numeric(y_1[1:i]),as.numeric(exog_var))
      
      
      for (j in (1:length(endog_exog))) {
        Z_1[i,(column+j-2)] <- endog_exog[j]
      }
      column <- column + i +3
    }
    
    
    x_1 <- exog_var
    
    x_i <- x_1
    
    
    Z_i <- matrix(0,T-2,n_inst+((T-2)*3))
    
    Z <- Z_1
    
    
    #We start at the column 30 and go by chunks of 29
    #but we will only loop over the Y chunks until 27 !
    #that is because yi27 will instrument delta_ei29
    for (sst in seq(31, nrow(data), T)) {
      
      #print(sst)
      #The second chunk for example goes from 30 to 
      #30 + 27 
      y_i <- data[sst:(sst+27),10]
      
      x_i <- data[(sst+2):((sst+2)+28),16:18]
      
      #The matrix of instruments has size 27 x 378 ! 
      #each time period we use more lags
      Z_i <- matrix(0,T-2,n_inst+((T-2)*3))
      
      column <-2 
      #Next time I will start in row 57 which is row 59 -2 
      #and so one the keeps increase by 1 for each loop, this -2 takes the name iter
      
      
      for (i in (1:(T-2))) {
        
        #Chunk is the bunch of yi that we are going to put in the Z_i matrix
        chunk <- as.numeric(y_1[1:i])
        #print(chunk)
        exog_var <- x_i[i,]
        
        
        endog_exog <- append(as.numeric(y_1[1:i]),as.numeric(exog_var))
        
        for (j in (1:length(endog_exog))) {
          Z_i[i,(column+j-2)] <- endog_exog[j]
        }
        column <- column + i +3
      }
      
      
      Z <- rbind(Z,Z_i)
      
    }
    
    
    
    n_inst_var <- N*(T-2)
    
    #Imagine our big H is very big -> we try to create 
    diagonal <- 2
    offdiagonal<- -1
    H <- matrix(0,n_inst_var,n_inst_var)
    diag(H) <- diagonal
    diag(H[-1,])<-offdiagonal
    diag(H[,-1])<-offdiagonal
    
    #this is a matrix with diagonal 2 and -1 on each side 
    
    data_reg <- na.omit(data)
    
    x1 <- as.matrix(na.omit(data_reg[,16:18]))
    x2 <- as.matrix(na.omit(data_reg[,19]))
    x<- cbind(x2,x1) #V1 like vendogemous
    
    data_reg <- na.omit(data)
    y <- as.matrix(na.omit(data_reg[,15]))
    
    W_notinv <- t(Z) %*% H %*% Z 
    
    #we have the big W optimal
    W_opt <- solve(W_notinv,tol = 1e-22)
    
    
    gamma <- solve(t(x) %*% Z %*% W_opt %*% t(Z) %*% x) %*%  t(x) %*% Z %*% W_opt %*% t(Z) %*% y
    
    y_pred <- y
    yhat   <- as.vector(x%*%gamma)
    res    <- y_pred-yhat
    
    n  <- length(y)
    k  <- ncol(x)
    df <- n-k
    
    sigma2 <- as.vector(t(res)%*%res)/df
    
    sigma2 <- sigma2/2
    
    var_gamma <- sigma2* solve(t(x) %*% Z %*% W_opt %*% t(Z) %*% x)
    
  
    
    stdvs_BGMM_sys  <-  sqrt(diag(var_gamma))
    stdvs_BGMM_sys
    
    tstats_BGMM_sys <- gamma/stdvs_BGMM_sys
    tstats_BGMM_sys <- c(tstats_BGMM_sys)
    
    n  <- length(y)
    k  <- ncol(x)
    df <- n-k
    
    pvals_GGMM_sys   <- 2*(1-pt(abs(tstats_BGMM_sys),df))
    
    
    #save output
    names(gamma) <- colnames(x)
    names(tstats_BGMM_sys) <- colnames(x)
    
    #creating the table
    coefs  <- round(gamma,3)
    stdvs  <- round(stdvs_BGMM_sys,3)
    tstats <- round(tstats_BGMM_sys,3)
    pvals  <- round(pvals_GGMM_sys,3)
    
    
    
    if(stack == 1){
      
      
      N <-   length(unique(data$state))
      
      
      #Let's loop over the states and create our big matrix
      
      
      yit <- seq(1,T-2)
      n_inst <- sum(yit)
      
      #The chunks have size 27 !!!
      #because yi27 is the last instrument of delta_ei29
      y_1 <- data[1:(T-2),10]
      
      
      
      #The number of columns of Z_1 is we would have without adding is 1 + 2 + .. + 27 = ninst
      #but now for each chunk 1 + 2 + ... + 27, we had the 3 explanatory variables as instrument
      
      Z_1 <- matrix(0,T-2,n_inst)
      column <-1 
      
      for (i in (1:(T-2))) {
        
        #Chunk is the bunch of yi that we are going to put in the Z_i matrix
        chunk <- as.numeric(y_1[1:i])
        
        chunk <- rev(chunk)
        
        column <- column + i -1
        
        for (j in (1:length(chunk))) {
          Z_1[i,column+j-1] <- chunk[j]
        }
        
      }
      
      Z <- Z_1
      Z_i <- matrix(0,T-2,n_inst+((T-2)*3))
      
      #We start at the column 30 and go by chunks of 29
      #but we will only loop over the Y chunks until 27 !
      #that is because yi27 will instrument delta_ei29
      for (sst in seq(31, nrow(data), T)) {
        
        #print(sst)
        #The second chunk for example goes from 30 to 
        #30 + 27 
        y_i <- data[sst:(sst+27),10]
        
        
        #The matrix of instruments has size 27 x 378 ! 
        #each time period we use more lags
        Z_i <- matrix(0,T-2,n_inst)
        column <-1 
        
        for (i in (1:(T-2))) {
          #print(i)
          chunk <- as.numeric(y_i[1:i])
          chunk <- rev(chunk)
          column <- column + i -1
          
          for (j in (1:length(chunk))) {
            Z_i[i,column+j-1] <- chunk[j]
          }
          
        }
        
        
        Z <- rbind(Z,Z_i)
        
      }
      
      
      #time to stack the matrix !
      
      Z[Z == 0] <- NA
      
      Z <- t(apply(Z,1,function(x){
        c(x[!is.na(x)],x[is.na(x)])}))
      
      
      Z <- Z[,colSums(is.na(Z))<nrow(Z)]
      
      
      Z[is.na(Z)] <- 0
      
    
      n_inst_var <- N*(T-2)
      
      #Imagine our big H is very big -> we try to create 
      diagonal <- 2
      offdiagonal<- -1
      H <- matrix(0,n_inst_var,n_inst_var)
      diag(H) <- diagonal
      diag(H[-1,])<-offdiagonal
      diag(H[,-1])<-offdiagonal
      
      #this is a matrix with diagonal 2 and -1 on each side 
      
      data_reg <- na.omit(data)
      
      x1 <- as.matrix(na.omit(data_reg[,16:18]))
      x2 <- as.matrix(na.omit(data_reg[,19]))
      x<- cbind(x1,x2) #V1 like vendogemous
      
      data_reg <- na.omit(data)
      y <- as.matrix(na.omit(data_reg[,15]))
      
      big_Z <- as.matrix(cbind(x1,Z))
      
      Z <- big_Z
      
      
      W_notinv <- t(Z) %*% H %*% Z 
      
      #we have the big W optimal
      W_opt <- solve(W_notinv,tol = 1e-22)
      
      
      gamma <- solve(t(x) %*% Z %*% W_opt %*% t(Z) %*% x) %*%  t(x) %*% Z %*% W_opt %*% t(Z) %*% y
      
      y_pred <- y
      yhat   <- as.vector(x%*%gamma)
      res    <- y_pred-yhat
      
      n  <- length(y)
      k  <- ncol(x)
      df <- n-k
      
      sigma2 <- as.vector(t(res)%*%res)/df
      
      sigma2 <- sigma2/2
      
      
      var_gamma <- sigma2* solve(t(x) %*% Z %*% W_opt %*% t(Z) %*% x)
      
      
      
      stdvs_BGMM_sys  <-  sqrt(diag(var_gamma))
      stdvs_BGMM_sys
      
      tstats_BGMM_sys <- gamma/stdvs_BGMM_sys
      tstats_BGMM_sys <- c(tstats_BGMM_sys)
      
      n  <- length(y)
      k  <- ncol(x)
      df <- n-k
      
      pvals_GGMM_sys   <- 2*(1-pt(abs(tstats_BGMM_sys),df))
      
      
      #save output
      names(gamma) <- colnames(x)
      names(tstats_BGMM_sys) <- colnames(x)
      
      #creating the table
      coefs  <- round(gamma,3)
      stdvs  <- round(stdvs_BGMM_sys,3)
      tstats <- round(tstats_BGMM_sys,3)
      pvals  <- round(pvals_GGMM_sys,3)
      
    
    }
  
    
    return(out_pvals_BGMM_sys = cbind(coefs,stdvs,tstats,pvals))
    
    
    
  } else {
    
    N <-   length(unique(data$state))
    T <-   length(unique(data$year)) #the real T
    
    
    #el is the number of maximum instruments per period
    el <- max_lag-1
    #Let's loop over the states and create our big matrix
    
    n_inst <- sum(seq(1,(max_lag-1))) + ((T-2)-length(seq(1,(max_lag-1))))*(max_lag-1)
    
    
    #The chunks have size 28 !!!
    #because yi28 is the last instrument of delta_ei30
    
    y_1 <- data[1:(T-2),10]
    
    #number of instruments and exog variables entering the z matrix
    n_var <- n_inst + (3*(T-2))
    
    Z_1 <- matrix(0,T-2,n_var)
    column <-2 
    
    for (i in (1:(T-2))) {
      
      #Chunk is the bunch of yi that we are going to put in the Z_i matrix
      chunk <- as.numeric(y_1[1:i])
      #print(chunk)
      exog_var <- data[i+2,16:18]
      
      if (length(chunk) < (max_lag)){
        
        endog_exog <- append(as.numeric(y_1[1:i]),as.numeric(exog_var))
        
        for (j in 1:length(endog_exog)){
          
          
          Z_1[i,(column+j-2)] <- endog_exog[j]
          
        }
        
        column <- column + i +3
        
      } else {
        
        el <- max_lag-1
        
        endog_exog <- append(as.numeric(tail(chunk,el)),as.numeric(exog_var))
        
        
        for (j in 1:length(endog_exog)){
          
          Z_1[i,column+j-2] <- endog_exog[j]
        }
        
        column <- column + el + 3
      }
    }
    
    
    Z <- Z_1
    
    #We start at the column 30 and go by chunks of 29
    #but we will only loop over the Y chunks until 27 !
    #that is because yi27 will instrument delta_ei29
    for (sst in seq(31, nrow(data), T)) {
      
      #print(sst)
      #The second chunk for example goes from 30 to 
      #30 + 27 
      y_i <- data[sst:(sst+27),10]
      
      x_i <- data[(sst+2):((sst+2)+28),16:18]
      
      #The matrix of instruments has size 27 x 378 ! 
      #each time period we use more lags
      Z_i <- matrix(0,T-2,n_inst+((T-2)*3))
      
      column <-2 
      #Next time I will start in row 57 which is row 59 -2 
      #and so one the keeps increase by 1 for each loop, this -2 takes the name iter
      
      
      for (i in (1:(T-2))) {
        
        
        #Chunk is the bunch of yi that we are going to put in the Z_i matrix
        chunk <- as.numeric(y_i[1:i])
        #print(chunk)
        exog_var <- x_i[i,]
        
        
        if (length(chunk) < (max_lag)){
          
          endog_exog <- append(as.numeric(y_i[1:i]),as.numeric(exog_var))
          
          for (j in 1:length(endog_exog)){
            
            Z_i[i,(column+j-2)] <- endog_exog[j]
            
          }
          
          column <- column + i +3
          
        } else {
          
          el <- max_lag-1
          
          endog_exog <- append(as.numeric(tail(chunk,el)),as.numeric(exog_var))
          
          
          for (j in 1:length(endog_exog)){
            
            Z_i[i,column+j-2] <- endog_exog[j]
          }
          
          column <- column + el + 3
        }
      }
      
      Z <- rbind(Z,Z_i)
      
    }
    
    
    #changing number of instrumeeents!
    n_inst_var <- N*(T-2)
    
    #Imagine our big H is very big -> we try to create 
    diagonal <- 2
    offdiagonal<- -1
    H <- matrix(0,n_inst_var,n_inst_var)
    diag(H) <- diagonal
    diag(H[-1,])<-offdiagonal
    diag(H[,-1])<-offdiagonal
    
    
    
    data_reg <- na.omit(data)
    
    x1 <- as.matrix(na.omit(data_reg[,16:18]))
    x2 <- as.matrix(na.omit(data_reg[,19]))
    x<- cbind(x2,x1) #V1 like vendogemous
    
    data_reg <- na.omit(data)
    y <- as.matrix(na.omit(data_reg[,15]))
    
    W_notinv <- t(Z) %*% H %*% Z 
    
    
    
    #we have the big W optimal
    W_opt <- solve(W_notinv,tol=1e-23)
    
    
    gamma <- solve(t(x) %*% Z %*% W_opt %*% t(Z) %*% x) %*%  t(x) %*% Z %*% W_opt %*% t(Z) %*% y
    
    y_pred <- y
    yhat   <- as.vector(x%*%gamma)
    res    <- y_pred-yhat
    
    n  <- length(y)
    k  <- ncol(x)
    df <- n-k
    
    sigma2 <- as.vector(t(res)%*%res)/df
    
    sigma2 <- sigma2/2
    
    var_gamma <- sigma2* solve(t(x) %*% Z %*% W_opt %*% t(Z) %*% x)
    
    
    stdvs_BGMM_sys  <-  sqrt(diag(var_gamma))
    stdvs_BGMM_sys
    
    tstats_BGMM_sys <- gamma/stdvs_BGMM_sys
    tstats_BGMM_sys <- c(tstats_BGMM_sys)
    
    
    pvals_GGMM_sys   <- 2*(1-pt(abs(tstats_BGMM_sys),df))
    
    
    #save output
    names(gamma) <- colnames(x)
    names(tstats_BGMM_sys) <- colnames(x)
    
    #creating the table
    coefs  <- round(gamma,3)
    stdvs  <- round(stdvs_BGMM_sys,3)
    tstats <- round(tstats_BGMM_sys,3)
    pvals  <- round(pvals_GGMM_sys,3)
    
    
    return(out_pvals_BGMM_sys = cbind(coefs,stdvs,tstats,pvals))
    
  }
  
  
  
}
