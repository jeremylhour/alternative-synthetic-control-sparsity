gamma <- function(beta,d,X){
  ### Calibration Loss function
  ### Jeremy L Hour
  ### 14 janvier 2016
  
  ### INPUTS :
  # beta : coefficients
  # d : treatment indicator
  # X : covariates
  
  X <- as.matrix(X)
  f <- (1-d)*exp(X%*%beta) - d * (X%*%beta)                                                                       
  return(mean(f,na.rm=T))
}