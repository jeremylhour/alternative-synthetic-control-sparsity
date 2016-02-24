msystem <- function(beta,theta,y,d,X){
  ### The m-system as in the paper
  
  w <- d - exp(X%*%beta)*(1-d)
  W <- diag(as.vector(w))
  Z <- cbind(y,X)
  
  m <- W %*% Z
  m[,1] <-  m[,1] - theta*d
  
  return(m)
}