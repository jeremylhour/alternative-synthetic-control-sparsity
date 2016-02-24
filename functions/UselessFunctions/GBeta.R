GBeta <- function(beta,y,d,X,n){
  ### The G_Beta matrix as in the paper
  
  n <- ncol(X)
  w <- as.vector(exp(X%*%beta)*(1-d))
  Z <- cbind(y,X)
  Z <- sweep(Z,MARGIN=1, w, '*')
  
  return(t(Z) %*% X /n)
}