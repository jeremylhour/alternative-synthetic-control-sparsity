VanillaOrthogonalityReg <- function(y,d,X,beta,method="WLSLasso",
                        c=1.1, nopenset=c(1), RescaleY=F,
                       maxIterLasso=1e6,tolLasso=1e-6,PostLasso=F,trace=F){
  ### Function to compute mu hat
  ### Jeremy L Hour
  ### 11 janvier 2016
  ### EDITED : 19 fevrier 2016
  
  ### INPUTS:
  # y : Outcome variable (not normalized)
  # d : Treatment indicator
  # X : Covariates (not normalized)
  # beta : estimate from the first step
  # method : c("WLSLasso","LinearOutcome")  
  # c : constant for the overall penalty level
  # nopenset : set of indices that should not be penalized
  #            Default is intercept not penalized
  # RescaleY : if TRUE rescale variable Y
  # maxIterPen : maximal number of iterations for penalty convergence
  # maxIterLasso : maximal number of iterations in Lasso procedure
  # tolLasso : tolerance for stopping criterion in Lasso minimization
  # PostLasso : if TRUE computes the PostLasso instead of Lasso
  # trace : if TRUE print convergence info
  
  # The first column of X must be the constant
  
  
  ### Load user-defined functions
  source("functions/LassoFISTA.R")
  
  ### Setting
  d <- as.matrix(d)
  y <- as.matrix(y)
  X <- as.matrix(X)
  
  n <- nrow(X)
  p <- ncol(X)
  
  ### Type of method to compute weights
  if(method == "WLSLasso"){ W=(1-d)*exp(X%*%beta) }
  if(method == "LinearOutcome"){ W=(1-d)*(sum(d)/sum(1-d)) }
  W <- as.vector(W)
  
  ### First step: Lasso
  
  # Overall penalty level
  g <- .1/log(max(p,n))
  lambda <- c*qnorm(1-.5*g/p)/sqrt(n)
  
  # Adjustment for rescaled program (only a numerical matter)
  if(RescaleY){
    y_tilde <- (y-mean(y))/sd(y)
    lambda_tilde <- lambda/sd(y)
  } else {
    y_tilde <- y
    lambda_tilde <- lambda
  }
  
  # Estimation parameters
  mu <- rep(0,p)
  
  # Compute Lasso estimate
  LassoEstim <- LassoFISTA(betaInit=mu,y_tilde,X,W=W,
                             nopen=nopenset,lambda=lambda_tilde,
                             tol=tolLasso,maxIter=maxIterLasso,trace=F)
  mu <- LassoEstim$beta
    
  # Obtain the estimates for the unscaled model
  if(RescaleY){
    muLasso <- sd(y)*mu
    muLasso[1] <- muLasso[1] + mean(y)
  } else {
    muLasso <- mu
  }
  
  # Get estimated active set
  SHat <- union(1,which(muLasso != 0)) 

  
  ### Second step: Post-Lasso
  muPL <- rep(0,p)
  if(PostLasso==T){
    OrthoReg <- lm(y ~ X[,SHat] - 1, weights=W)
    muPL[SHat] <- OrthoReg$coefficients
    muPL[is.na(muPL)] <- 0
  }
  
  # Return objects
  return(list(method=method,
              lambda=lambda,
              SHat=SHat,
              muPL=muPL,
              muLasso=c(muLasso)
  ))
  
  
}
