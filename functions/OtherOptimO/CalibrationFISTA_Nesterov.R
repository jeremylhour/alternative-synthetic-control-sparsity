CalibrationFISTA_Nesterov <- function(betaInit=rep(0,ncol(X)),d,X,lambda,
                                     tol=1e-6,maxIter=1000,maxLineSearch=10,trace=F){
  ### Computes Lasso solution for Calibration objective function
  ### with Nesterov second method and line search
  ### Jeremy L Hour
  ### 14 janvier 2016
  
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
  
  t_start <- Sys.time()
  ### Set Algo. Values
  t <- .1
  b <- .5
  beta <- betaInit
  v <- beta
  
  f_path <- NULL
  
  k <- 0
  
  repeat{
    k <- k+1
    theta <- 2/(k+1)
    betaO <- beta
    
    z <- (1-theta)*beta + theta*v
    
    j <- 0
    repeat{
      j <- j+1
      v <- prox(v - t * gammagrad(z,d,X)/theta, lambda*t/theta)
      beta <- (1-theta)*beta + theta*v
      if(gamma(beta,d,X) < gamma(z,d,X) + t(gammagrad(z,d,X))%*%(beta-z) + sum((beta-z)^2)/(2*t) || j > maxLineSearch) break
      t <- b*t
    }
    
    f_path <- c(f_path,LassoObj(beta,d,X,lambda))
    
    if(trace==T & k%%5==0){ print(paste("Objective Func. Value at step",k,":",LassoObj(beta,d,X,lambda))) }
    
    if(abs( LassoObj(beta,d,X,lambda) - LassoObj(betaO,d,X,lambda) ) < tol || k > maxIter) break
  }
  
  
  if(k > maxIter){ print("Lasso minimization did not converge.")
                   convergence=-999
  } else {
    convergence=0
  }
  
  return(list(beta=beta,
              lambda=lambda,
              value=LassoObj(beta,d,X,lambda),
              nbIter=k,
              convergence=convergence,
              duration=Sys.time()-t_start,
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