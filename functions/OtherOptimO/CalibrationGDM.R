CalibrationGDM <- function(betaInit=rep(0,ncol(X)),d,X,
                             tol=1e-6,maxIter=1000,maxLineSearch=10,trace=F){
  ### Computes solution for Calibration objective function without penalty
  ### with Gradient Descent method
  ### Jeremy L Hour
  ### 15 janvier 2016
  
  ## INPUTS:
  # betaInit : Starting value for the coefficient
  # y : response
  # X : Covariates
  # W : Weights for each observation in the Least Square Objective
  # lambda : Lasso penalty parameter
  
  # tol : acceptable difference between the value of objective function at two iterations
  # maxIter : maximal number of iterations
  # maxLineSearch : maximal number of iterations for line search step
  # trace : if TRUE print convergence info
  
  ### Load user-defined functions
  source("functions/gamma.R")
  source("functions/gammagrad.R")
  
  X <- as.matrix(X)
  
  ### Set Algo. Values
  beta <- betaInit
  t <- 1
  b <- .99
  a <- .33
  
  k <- 0
  
  repeat{
    k <- k+1
    delta <- - gammagrad(beta,d,X)
    
    # Line Search
    j <- 0
    repeat{
      j <- j+1
      t <- b*t
      if(gamma(beta+t*delta,d,X) < gamma(beta,d,X) + a*t*t(gammagrad(beta,d,X))%*%delta || j > maxLineSearch) break
    }
    
    # Update
    betaO <- beta
    beta <- betaO + t*delta
    
    if(trace==T & k%%5==0) print(paste("Objective Func. Value at step",k,":",gamma(beta,d,X)))
    
    if(abs(gamma(beta,d,X)-gamma(betaO,d,X)) < tol || k > maxIter) break 
  }
  
  if(k > maxIter) print("Reach max. number of iterations reach in Lasso minimization.")
  
  return(list(beta=beta,
              value=gamma(beta,d,X),
              nbIter=k))
}