OrthogonalityReg <- function(y,d,X,beta,method="LinearOutcome",c=1.1,
                             maxIterPen=100,PostLasso=F,trace=F){
  ### Function to compute beta_0, the covariate balancing weights
  ### Jeremy L Hour
  ### 11 janvier 2016
  
  ### INPUTS:
  # y : Outcome variable
  # d : Treatment indicator
  # X : Covariates
  # beta : estimate from the first step
  # method : c("CHS2015","WLSLasso","LinearOutcome")  
  # c : constant for the overall penalty level
  # maxIterPen : maximal number of iterations for penalty convergence
  # PostLasso : if TRUE computes the PostLasso instead of Lasso
  # trace : if TRUE print convergence info
  
  # The first column of X must be the constant

  
  ### Load packages
  library("lbfgs")
  
  ### Load user-defined functions
  source("functions/LeastSqgrad.R")
  source("functions/LeastSq.R")
  
  ### Setting
  d <- as.matrix(d)
  y <- as.matrix(y)
  X <- as.matrix(X)
  
  n <- sum(1-d)
  p <- ncol(X)
  
  ### Type of method
  if(method %in% c("WLSLasso","CHS2015")){ W=(1-d) * exp(X%*%beta)/sum(d) }
  if(method == "LinearOutcome"){ W=1-d }
  
  ### For CHS2015 verify a few things:
  if(method == "CHS2015" & p > n){ stop("Not enough data to use CHS 2015 orthogonalisation.")}
  
  ### First step: Lasso
  
  # Overall penalty level
  g <- .1/log(max(p,n))
  lambda <- c*qnorm(1-.5*g/p)/sqrt(n)
  
  # Penalty loadings: get a preliminary estimate
  Psi <- diag(as.vector(sqrt( t((y-mean(y))^2) %*% X^2 / n )))
  
  # Estimation parameters
  v <- .01 # Stopping rule
  k <- 0
  mu <- rep(0,p)
  
  if(method!="CHS2015"){
    
    # Lasso estimate
    repeat{
      k <- k+1
      
      # Compute Lasso estimate
      LassoEstim <- lbfgs(LeastSq, LeastSqgrad, Psi%*%mu,y=y,X=X%*%solve(Psi),W=W,
                          orthantwise_c=lambda,orthantwise_start=1,
                          invisible=1)
      mu <- solve(Psi) %*% LassoEstim$par
      
      # Update penalty loadings
      PrePsi <- Psi
      Psi <- diag(as.vector(sqrt( t((y-X%*%mu)^2) %*% X^2 / n )))
      
      if(trace & k%%5==0){
        print(paste("Max. pen. loading diff at Lasso Iteration nb.",k,":",max(abs(diag(Psi-PrePsi))))) 
      }
      
      if(k > maxIterPen || max(abs(diag(Psi-PrePsi))) < v) break
    }
    
    if(k > maxIterPen & method!="CHS2015") print("Reach max. number of iterations in penalty estimation.")
    
    muLasso <- mu 
    # Get estimated active set
    SHat <- union(1,which(muLasso != 0)) # Always put a constant
    
  }
  
  
  ### Adjustement for CHS 2015 method
  if(method == "CHS2015"){SHat = 1:p}
  
  
  ### Second step: Post-Lasso
  if(PostLasso==T || method == "CHS2015"){
    OrthoReg <- lm(y ~ X[,SHat] - 1, weights=W)
    muPL <- rep(0,p)
    muPL[SHat] <- OrthoReg$coefficients
    muPL[is.na(muPL)] <- 0
  }
  
  
  if(k > maxIterPen){ print("Penalty estimation did not converge.")
                      cvg=-999
  } else {
    cvg=0
  }
  
  if(PostLasso==T || method == "CHS2015"){
    return(list(method=method,SHat=SHat,muPL=muPL,muLasso=c(muLasso),lambda=lambda,nbIter=k,convergence=cvg))
  } else {
    return(list(method=method,SHat=SHat,muLasso=c(muLasso),lambda=lambda,nbIter=k,convergence=cvg))
  }
  
}