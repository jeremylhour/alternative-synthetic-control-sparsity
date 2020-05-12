#' Function to compute beta, the covariate balancing parameter
#' 
#' First step of the BEAST estimator. Uses the lbfgs function to perform L1-penalised minimization.
#' A constant must be included as the first column in X and is not penalized.
#' Last edited: 12 mai 2020.
#' 
#' @param d Treatment indicator.
#' @param X Matrix of covariates.
#' @param c Constant for the overall penalty level.
#' @param maxIterPen Maximal number of iterations for penalty estimation.
#' @param trace if TRUE print convergence info.
#' @param PostLasso if TRUE computes the Post-Lasso solution.
#' @param maxIter Maximal number of iteration for Lasso program.
#' 
#' @return SHat Set of indices of non-zero elements in estimated beta.
#' @return betaLasso Lasso solution.
#' @return betaPL Post-Lasso solution.
#' @return lambda Overall penalty level.
#' @return psi Covariate-specific penalty loadings.
#' @return nbIter Number of iterations for penalty level estimation.
#' @return convergence 0 if convergence, -999 if not.
#' 
#' @author Jeremy L'Hour


CalibrationLasso <- function(d,X,c=1.1,maxIterPen=100,PostLasso=F,trace=F){
  ### Setting
  d <- as.matrix(d)
  X <- as.matrix(X)
  
  n <- nrow(X)
  p <- ncol(X)
  
  ### First step: Lasso
  
  # Overall penalty level
  g <- 10/sqrt(log(max(p,n)))
  lambda <- c*qnorm(1-.5*g/p)/sqrt(n)
  
  # Penalty loadings: get a preliminary estimator
  beta <- c(log(sum(d)/(sum(1-d))),rep(0,p-1))
  Psi <- diag(as.vector(sqrt( t((d - exp(X%*%beta)*(1-d))^2) %*% X^2 / n )))
  
  # Estimation parameters
  v <- .01 # Stopping rule
  k <- 0
  
  # Lasso estimate
  repeat{
    k <- k+1
    LassoEstim <- lbfgs(gamma, gammagrad, Psi%*%beta, d=d, X=X%*%solve(Psi),
                        orthantwise_c=lambda,orthantwise_start=1,
                        invisible=1)
    beta <- solve(Psi) %*% LassoEstim$par
    
    # Update penalty loadings
    PrePsi <- Psi
    Psi <- diag(as.vector(sqrt( t((d - exp(X%*%beta)*(1-d))^2) %*% X^2 / n )))
    
    if(trace & k%%5==0) print(paste("Max. pen. loading diff at Lasso Iteration nb.",k,":",max(abs(diag(Psi-PrePsi)))))
    if(any(exp(X%*%beta)==Inf)) break
    if(k > maxIterPen || max(abs(diag(Psi-PrePsi))) < v) break
  }
  
  betaLasso <- beta 
  SHat <- union(1,which(betaLasso != 0)) # Always put a constant
  
  
  ### Second step: Post-Lasso
  betaPL <- rep(0,p)
  if(PostLasso==T){
    PL <- lbfgs(gamma, gammagrad, beta[SHat], d=d, X=X[,SHat],
                orthantwise_c=0, invisible=1)
    betaPL[SHat] <- PL$par
  }
  
  
  if(k > maxIterPen || any(exp(X%*%beta)==Inf)){
    cvg=-999
    if(trace) print("Penalty estimation did not converge.")
  } else {
    cvg=0
  }
  
  return(list(SHat=SHat,
              betaPL=betaPL,
              betaLasso=c(betaLasso),
              lambda=lambda, 
              nbIter=k, 
              convergence=cvg,
              psi=diag(Psi)))
  
  
}

##################################
##################################
##################################
### Define Auxiliary Functions ###
##################################
##################################
##################################

gamma <- function(beta,d,X){
  X <- as.matrix(X)
  f <- (1-d)*exp(X%*%beta) - d * (X%*%beta)                                                                       
  return(mean(f,na.rm=T))
}

gammagrad <- function(beta,d,X){
  X <- as.matrix(X)
  g <- as.vector(( (1-d)*exp(X%*%beta) - d )) * X
  return(as.vector(apply(g,2,mean,na.rm=T)))
}