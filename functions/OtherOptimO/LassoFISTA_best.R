LassoFISTA <- function(betaInit=rep(0,ncol(X)),y,X,W=rep(1,nrow(X)),lambda,
                       tol=1e-6,maxIter=1000,trace=F){
  ### Computes Lasso solution for Least Square objective function
  ### with FISTA method
  ### Jeremy L Hour
  ### 15 janvier 2016
  
  ## INPUTS:
  # betaInit : Starting value for the coefficient
  # y : response
  # X : Covariates
  # W : Weights for each observation in the Least Square Objective
  # lambda : Lasso penalty parameter
  # tol : acceptable difference between the value of objective function at two iterations
  # K : maximal number of iterations
  # trace : if TRUE print convergence info
  
  # The first column of X should be the constant
  
  ### Load user-defined functions
  source("functions/LeastSqgrad.R")
  source("functions/LeastSq.R")
  
  ### Set Algo. Values
  eta <- 1/max(2*eigen(t(X)%*%X)$values/nrow(X))
  theta <- 1
  thetaO <- theta
  beta <- betaInit
  v <- beta
  
  k <- 0
  repeat{
    k <- k+1
    
    thetaO <- theta
    theta <- (1+sqrt(1+4*thetaO^2))/2
    delta <- (1-thetaO)/theta
    
    betaO <- beta
    beta <- prox(v - eta*LeastSqgrad(v,y,X,W), lambda*eta)
    
    v <- (1-delta)*beta + delta*betaO
    
    if(trace==T & k%%5 == 0){ print(paste("Objective Func. Value at step",k,":",LassoObj(beta,y,X,W,lambda))) }
    
    if(abs( LassoObj(beta,y,X,W,lambda) - LassoObj(betaO,y,X,W,lambda) ) < tol || k > maxIter) break

  }
  
  if(k > maxIter) print("Reach max. number of iterations reach in Lasso minimization.")

  return(list(beta=beta,
              lambda=lambda,
              value=LassoObj(beta,y,X,W,lambda),
              nbIter=k))
}


### Define auxiliary functions
prox <- function(x,lambda){
  y <- (abs(x)-lambda)*(abs(x)-lambda > 0) * sign(x)
  y[1] <- x[1] # Do not penalize constant
  return(y)
}

LassoObj <- function(beta,y,X,W,lambda){
  LeastSq(beta,y,X,W) + lambda*sum(abs(beta))
}