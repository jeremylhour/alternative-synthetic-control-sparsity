CalibrationNesterov <- function(betaInit,d,X,lambda,tol=1e-6,maxIter=1000,trace=F){
  ### Computes Lasso Calibration solution with Nesterov method
  ### Jeremy L Hour
  ### 13 janvier 2016
  
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
  
  source("functions/gammagrad.R")
  source("functions/gamma.R")
  
  
  t <- .001 # step size
  beta <- betaInit
  v <- beta
  
  f_path <- vector(length=maxIter)
  
  k <- 0
  
  repeat{
    k <- k+1
    
    theta <- 2/(k+1)
    betaO <- beta
    u <- (1-theta)*beta + theta * v
    v <- prox(v - (t/theta) * gammagrad(u,d,X), lambda*t/theta)
    beta <- (1-theta)*beta + theta * v
    
    if(trace==T){ print(paste("Objective Func. Value at step",k,":",LassoObj(beta,d,X,lambda))) }
    f_path[k] <- LassoObj(beta,d,X,lambda)
    if(abs(gamma(beta,d,X)-gamma(betaO,d,X) + lambda*sum(abs(beta)) - lambda*sum(abs(betaO))) < tol || k > maxIter) break 
  }
  
  
  return(list(beta=beta,
              lambda=lambda,
              value=gamma(beta,d,X),
              f_path=f_path))
}

### Define auxiliary functions
prox <- function(x,lambda){
  (abs(x)-lambda)*(abs(x)-lambda > 0) * sign(x)
}

LassoObj <- function(beta,d,X,lambda){
  gamma(beta,d,X) + lambda*sum(abs(beta))
}