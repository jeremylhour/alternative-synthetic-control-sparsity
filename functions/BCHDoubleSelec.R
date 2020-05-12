#' BCH (2014) Double-selection procedure
#' 
#' Performs the double selection procedure of Belloni, Chernozhukov and Hansen (JEP, 2014)
#'  A constant must be included as the first column in X.
#'  Last edited: 24 fevrier 2016.
#' 
#' @param y Outcome variable, not normalized.
#' @param X Matrix of covariates.  
#' @param cd Constant for the overall penalty level in Treatment Selection.
#' @param cy Constant for the overall penalty level in Outcome Selection.
#' @param nopenset Set of indices that should not be penalized. Default is intercept penalized.
#' @param RescaleY if TRUE rescale variable y.
#' @param maxIterPen Maximal number of iterations for penalty estimation.
#' @param maxIterLasso Maximal number of iterations in Lasso procedure.
#' @param tolLasso Tolerance for stopping criterion in Lasso minimization.
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



BCHDoubleSelec <- function(y,d,X,cd=1.1,cy=1.1,
                             nopenset=c(1),RescaleY=F,
                             maxIterPen=1e4,maxIterLasso=1e6,tolLasso=1e-6,trace=F){
  
  ### Setting
  d <- as.matrix(d)
  y <- as.matrix(y)
  X <- as.matrix(X)
  n <- nrow(X)
  p <- ncol(X)
  
  ### First step: Selection on the Outcome Equation
  
  # Overall penalty level
  g <- 10/sqrt(log(max(p,n)))
  lambda <- cy*2*qnorm(1-.5*g/p)/sqrt(n)
  
  # Adjustment for rescaled program (only a numerical matter)
  if(RescaleY){
    y_tilde <- (y-mean(y))/sd(y)
    lambda_tilde <- lambda/sd(y)
  } else {
    y_tilde <- y
    lambda_tilde <- lambda
  }
  
  # Penalty loadings: get a preliminary estimate
  m_y <- mean(y_tilde)
  Psi <- diag(as.vector(sqrt( t((y_tilde-m_y)^2) %*% X^2 / n )))
  
  # Estimation parameters
  v <- .01 # Stopping rule
  k <- 0
  muy <- rep(0,p)
  
  repeat{
    k <- k+1
    
    # Compute Lasso estimate
    LassoEstim <- LassoFISTA(betaInit=Psi%*%muy,y_tilde,X%*%solve(Psi),W=rep(1,n),
                             nopen=nopenset,lambda=lambda_tilde,
                             tol=tolLasso,maxIter=maxIterLasso,trace=F)
    muy <- solve(Psi)%*%LassoEstim$beta
    SHaty <- union(1,which(muy != 0)) 
    PostLasso <- lm(y_tilde ~ X[,SHaty]-1)
    e <- PostLasso$residuals
    
    # Update penalty loadings
    PrePsi <- Psi
    Psi <- diag(as.vector(sqrt( t(e^2) %*% X^2 / n )))
    
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
  
  if(k > maxIterPen){
    cvg=-999
    print("Penalty estimation did not converge.")
  }
  
  if(trace) print("Outcome selection step over.")
  
  
  ### Second step: Selection on the Treatment Equation
  
  # Overall penalty level
  lambda <- cd*2*qnorm(1-.5*g/p)/sqrt(n)
  
  # Penalty loadings: get a preliminary estimate
  m_d <- mean(d)
  Psi <- diag(as.vector(sqrt( t((d-m_d)^2) %*% X^2 / n )))
  
  # Estimation parameters
  v <- .01 # Stopping rule
  k <- 0
  mud <- rep(0,p)
  
  repeat{
    k <- k+1
    
    # Compute Lasso estimate
    LassoEstim <- LassoFISTA(betaInit=Psi%*%mud,d,X%*%solve(Psi),W=rep(1,n),
                             nopen=nopenset,lambda=lambda,
                             tol=tolLasso,maxIter=maxIterLasso,trace=F)
    mud <- solve(Psi)%*%LassoEstim$beta
    SHatd <- union(1,which(mud != 0)) 
    PostLasso <- lm(d ~ X[,SHatd]-1)
    e <- PostLasso$residuals
    
    # Update penalty loadings
    PrePsi <- Psi
    Psi <- diag(as.vector(sqrt( t(e^2) %*% X^2 / n )))
    
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
  
  if(k > maxIterPen){
    cvg=-999
    print("Penalty estimation did not converge.")
  }
  
  if(trace) print("Treatment selection step over.")
  
  
  ### Third step step: Post-Lasso
  DBSelecSet <- union(SHaty,SHatd)
  DBPostSelec <- lm(y ~ d + X[,DBSelecSet]-1)
  
  ### Fourth step : Compute the Variance
  PolicyReg <- lm(d ~ X[,DBSelecSet])
  v <- PolicyReg$residuals
  xi <- sqrt(n/(n-length(DBSelecSet)-1))*DBPostSelec$residuals
  
  sigma2 <- mean(v^2*xi^2) / mean(v^2)^2
  sigma <- sqrt(sigma2/n)
  
  # Return objects
  return(list(SHaty=SHaty,
              SHatd=SHatd,
              theta=coef(DBPostSelec)["d"],
              sigma=sigma,
              convergence=cvg
  )) 
}
