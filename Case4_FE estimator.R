####################################################################
##################              FE estimator             ##################
####################################################################

FE = function (y,x,T,TFE) 
###build matrix D
{
  N   <-nT/T
  D   <- matrix(0, NT, N)
  M_D <- diag(NT)-D %*% slove(t(D) %*% D) %*% t(D) ###Residual maker matrix
  ##with estimator follow slides p21 equation
  b   <- slove(t(x) %*% M_D %*% x) %*% t(x) %*% M_D  %*% y
  ##constant alpha
  a   <- slove(t(D) %*% D) %*% t(D) %*% (y-x %*% b)
  coefs_bw_fe<-b
  coefs_fe<-a
  ##FE variance-covariance matrix
  df    <- nT-N-K
  res_fe<- M_D %*% y - x %*% b
  sigma_ep_2 <- (t(res_fe) %*% res_fe)/ df
  vcm_bw <- sigma_ep_2 %*% slove(t(x) %*% M_D %*% x)
  std_fe<-sqrt(vcm_bw)
  tstats_fe = coefs_bw_fe/std_fe
  pvals_fe = 2*(1-pt(abs(tstats_fe),df))

}