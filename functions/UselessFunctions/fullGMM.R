fullGMM <- function(par,Z){
  ### Gradient of the loss function 
  ### Z is such that the first column is the response, then the outcome and X the covariates
  ### Similar to the m-system
  
  d <- as.vector(Z[,1])
  y <- as.vector(Z[,2])
  X <- as.matrix(Z[,3:ncol(Z)])
  
  w <- d - exp(X%*%matrix(par[2:length(par)], ncol=1))*(1-d)
  W <- diag(as.vector(w))
  Z <- cbind(y,X)
  
  m <- W %*% Z
  m[,1] <-  m[,1] - par[1]*d
  
  return(m)
}