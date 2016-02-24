ImmunizedATTVariance <- function(y,d,X,beta,mu,Immunity=T){
  ### Function to compute asymptotic variance 
  ### of the Immunized estimator
  ### Jeremy L Hour
  ### 5 fevrier 2016
  
  ### INPUTS:
  # y : Outcome variable
  # d : Treatment indicator
  # X : Covariates
  # beta : estimate from the first step
  # mu : Orthogonality parameter
  # Immunity : should be TRUE if Immunized and FALSE if Naive plug-in
  
  
  ### Setting
  d <- as.matrix(d)
  y <- as.matrix(y)
  X <- as.matrix(X)
  
  eps <- y
  if(Immunity==T){ eps <- y - X%*%mu }
  
  pi <- mean(d)
  theta <- mean((d - (1-d)*exp(X%*%beta)) * eps)  / pi
  
  psi <- (d - (1-d)*exp(X%*%beta)) * eps - d*theta
  
  return(mean(psi^2)/pi^2)
}