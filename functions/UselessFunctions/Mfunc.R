Mfunc <- function(par,data,beta){
  
  theta <- par[1]
  mu <- par[-1]
  
  y <- data$y
  d <- data$d
  X <- data$X
  
  m <- msystem(beta,theta,y,d,X)
  
  sqM <- mean(m %*% mu)^2
  
  return(sqM)
}