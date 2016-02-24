LeastSqgrad <- function(mu,y,X,W){
  ### Gradient of function to be minimized
  ### Jeremy L Hour
  ### 11 janvier 2016
  
  ### INPUTS:
  # y : Outcome variable
  # X : Covariates
  # W : (1-d)*exp(X %*% beta)

  X <- as.matrix(X)
  n <- nrow(X)
  P <- -2*W*(y - X%*%mu)
  return(
    as.vector(t(P) %*% X / n)
    )
}