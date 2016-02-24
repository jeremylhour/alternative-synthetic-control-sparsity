CalibrationFISTA <- function(betaInit,d,X,lambda,tol=1e-6,maxIter=1000,trace=F){
  ### Computes Lasso solution for Calibration objective function
  ### with FISTA method with line search (method 1)
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
  source("functions/gamma.R")
  source("functions/gammagrad.R")
  
  X <- as.matrix(X)
  
  ### Set Algo. Values
  L <- 2 * max(eigen(t(X)%*%X)$values) / nrow(X)
  t <- 1/L
  b <- .99
  beta <- betaInit
  v <- beta
  
  f_path <- vector(length=maxIter)
  
  k <- 0
  
  repeat{
    k <- k+1
    theta <- 2/(k+1)
    betaO <- beta
    
    u <- (1-theta)*beta + theta*v
    beta <- prox(u - t * gammagrad(u,d,X), lambda*t)
    
    j <- 0
    repeat{
      j <- j+1
      t <- b*t
      beta <- prox(u - t * gammagrad(u,d,X), lambda*t)
      if(gamma(beta,d,X) < gamma(u,d,X) + t(gammagrad(u,d,X))%*%(beta-u) + sum((beta-u)^2)/(2*t) || j > 9){
        t <- 1/L
        break
      } 
    }
    
    v <- betaO + (beta - betaO)/theta
    
    if(trace==T){ print(paste("Objective Func. Value at step",k,":",LassoObj(beta,d,X,lambda))) }
    
    f_path[k] <- LassoObj(beta,d,X,lambda)
    
    if(abs( LassoObj(beta,d,X,lambda) - LassoObj(betaO,d,X,lambda) ) < tol || k > maxIter) break 
  }
  
  return(list(beta=beta,
              lambda=lambda,
              value=LassoObj(beta,d,X,lambda),
              f_path=f_path))
}



### AUXILIARY FUNCTIONS

# Proximal operator for lambda * l1norm(coeff)
prox <- function(x,lambda){
  (abs(x)-lambda)*(abs(x)-lambda > 0) * sign(x)
}

# Lasso Objective function
LassoObj <- function(beta,d,X,lambda){
  gamma(beta,d,X) + lambda*sum(abs(beta))
}