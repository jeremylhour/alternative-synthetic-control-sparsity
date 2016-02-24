gammagrad <- function(beta,d,X){
  ### Calibration Loss function gradient
  ### Jeremy L Hour
  ### 14 janvier 2016
  
  ### INPUTS :
  # beta : coefficients
  # d : treatment indicator
  # X : covariates
  
  X <- as.matrix(X)
  g <- as.vector(( (1-d)*exp(X%*%beta) - d )) * X
  
  return(as.vector(apply(g,2,mean,na.rm=T)))
}