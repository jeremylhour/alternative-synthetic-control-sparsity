LeastSq <- function(mu,y,X,W){
  ### Loss function to be minimized
  ### Jeremy L Hour
  ### 11 janvier 2016
  
  ### INPUTS:
  # y : Outcome variable
  # X : Covariates
  # c : constant for the overall penalty level
  # K : maximal number of iterations
  # trace : if TRUE print convergence info
  
  X <- as.matrix(X)
  f <- W*(y - X%*%mu)^2
  return(mean(f))
}