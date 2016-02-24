LassoNesterov <- function(betaInit,y,X,W=rep(1,nrow(X)),lambda,tol=1e-6,maxIter=1000,trace=F){
  ### Computes Lasso solution for Least Square objective function
  ### with Nesterov second's method
  ### Jeremy L Hour
  ### 14 janvier 2016
  
  ## INPUTS:
  # betaInit : Starting value for the coefficient
  # y : response
  # X : Covariates
  # W : Weights for each observation in the Least Square Objective
  # lambda : Lasso penalty parameter
  # tol : acceptable difference between the value of objective function at two iterations
  # K : maximal number of iterations
  # trace : if TRUE print convergence info
  
  ### Load user-defined functions
  source("functions/LeastSqgrad.R")
  source("functions/LeastSq.R")
  
  ### Set Algo. Values
  L <- 2 * max(eigen(t(X)%*%X)$values) / nrow(X)
  t <- 1/L
  beta <- betaInit
  v <- beta
  
  f_path <- vector(length=maxIter)
  
  k <- 0
  
  repeat{
    k <- k+1
    theta <- 2/(k+1)
    betaO <- beta
    u <- (1-theta)*beta + theta*v
    v <- prox(v - (t/theta) * LeastSqgrad(u,y,X,W), lambda*t/theta)
    beta <- (1-theta)*beta + theta*v
    
    if(trace==T){ print(paste("Objective Func. Value at step",k,":",LassoObj(beta,y,X,W,lambda))) }
    
    f_path[k] <- LassoObj(beta,y,X,W,lambda)
    
    if(abs( LassoObj(beta,y,X,W,lambda) - LassoObj(betaO,y,X,W,lambda) ) < tol || k > maxIter) break 
  }
  
  return(list(beta=beta,
              lambda=lambda,
              value=LeastSq(beta,y,X,W),
              f_path=f_path))
}


### Define auxiliary functions
prox <- function(x,lambda){
  (abs(x)-lambda)*(abs(x)-lambda > 0) * sign(x)
}

LassoObj <- function(beta,y,X,W,lambda){
  LeastSq(beta,y,X,W) + lambda*sum(abs(beta))
}