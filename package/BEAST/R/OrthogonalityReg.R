#' Function to compute mu, the covariate balancing weights
#' 
#' Second step of the BEAST estimator. Uses the LassoFISTA function to perform L1-penalised minimization.
#'  A constant must be included as the first column in X.
#'  Last edited: 18 avril 2016.
#' 
#' @param y Outcome variable, not normalized.
#' @param X Matrix of covariates.
#' @param beta Calibrating parameter estimate from the first step.
#' @param method One of "WLSLasso" or "LinearOutcome".  
#' @param c Constant for the overall penalty level.
#' @param nopenset Set of indices that should not be penalized. Default is intercept penalized.
#' @param RescaleY if TRUE rescale variable y.
#' @param maxIterPen Maximal number of iterations for penalty estimation.
#' @param maxIterLasso Maximal number of iterations in Lasso procedure.
#' @param tolLasso Tolerance for stopping criterion in Lasso minimization.
#' @param PostLasso if TRUE computes the PostLasso solution.
#' @param trace if TRUE print convergence info.
#' 
#' @return SHat Set of indices of non-zero elements in estimated mu.
#' @return muLasso Lasso solution.
#' @return muPL Post-Lasso solution.
#' @return lambda Overall penalty level.
#' @return psi Covariate-specific penalty loadings.
#' @return nbIter Number of iterations for penalty level estimation.
#' @return convergence 0 if convergence, -555 if not because of Lasso minimization, -999 if not because of penalty estimation.
#' 
#' @seealso \code{\link{LassoFISTA}}
#' 
#' @author Jeremy Lhour



OrthogonalityReg <- function(y,d,X,beta,method="WLSLasso",
                             c=1.1,nopenset=c(1),RescaleY=F,
                             maxIterPen=1e4,maxIterLasso=1e6,tolLasso=1e-6,PostLasso=F,trace=F){
  
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
  
  # Penalty loadings: get a preliminary estimate
  m_y <- c(t(W)%*%y_tilde/sum(W))
  Psi <- diag(as.vector(sqrt( t(W*(y_tilde-m_y)^2) %*% (diag(sqrt(W))%*%X)^2 / n )))
  
  # Estimation parameters
  v <- .01 # Stopping rule
  k <- 0
  mu <- rep(0,p)
  
  repeat{
    k <- k+1
    
    # Compute Lasso estimate
    LassoEstim <- LassoFISTA(betaInit=Psi%*%mu,y_tilde,X%*%solve(Psi),W=W,
                             nopen=nopenset,lambda=lambda_tilde,
                             tol=tolLasso,maxIter=maxIterLasso,trace=F)
    mu <- solve(Psi)%*%LassoEstim$beta
    
    # Update penalty loadings
    PrePsi <- Psi
    Psi <- diag(as.vector(sqrt( t(W*(y_tilde-X%*%mu)^2) %*% (diag(sqrt(W))%*%X)^2 / n )))
    
    # Trace showing
    if(trace & k%%5==0){
      print(paste("Max. pen. loading diff at Lasso Iteration nb.",k,":",max(abs(diag(Psi-PrePsi))))) 
    }
    
    # Stopping rules
    if(k > maxIterPen || max(abs(diag(Psi-PrePsi))) < v || LassoEstim$convergenceFISTA==-555){
      cvg <- LassoEstim$convergenceFISTA
      break
    } 
  }
  
  cvg = 0
  if(k > maxIterPen){
    cvg=-999
    if(trace) print("Penalty estimation did not converge.")
  }
  
  # Obtain the estimates for the unscaled model
  if(RescaleY){
    muLasso <- sd(y)*mu
    muLasso[1] <- muLasso[1] + mean(y)
  } else {
    muLasso <- mu
  }
  
  # Get estimated active set
  if(cvg==0){
    SHat <- union(1,which(muLasso != 0)) 
  } else {
    SHat <- NA
  }
  
  ### Second step: Post-Lasso
  muPL <- rep(0,p)
  if(PostLasso==T & cvg==0){
    OrthoReg <- lm(y ~ X[,SHat] - 1, weights=W)
    muPL[SHat] <- OrthoReg$coefficients
    muPL[is.na(muPL)] <- 0
  }
  
  # Return objects
  return(list(lambda=lambda,
              SHat=SHat,
              muPL=muPL,
              muLasso=c(muLasso),
              nbIter=k,
              convergence=cvg,
              psi=diag(Psi)
  )) 
}
