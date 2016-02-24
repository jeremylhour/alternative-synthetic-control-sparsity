OrthogonalityReg <- function(y,d,X,beta,method="WLSLasso",c=1.1,
                             maxIterPen=100,maxIterLasso=1e6,tolLasso=1e-6,PostLasso=F,trace=F){
  ### Function to compute beta_0, the covariate balancing weights
  ### Jeremy L Hour
  ### 11 janvier 2016
  
  ### INPUTS:
  # y : Outcome variable
  # d : Treatment indicator
  # X : Covariates
  # beta : estimate from the first step
  # method : c("WLSLasso","LinearOutcome")  
  # c : constant for the overall penalty level
  # maxIterPen : maximal number of iterations for penalty convergence
  # maxIterLasso : maximal number of iterations in Lasso procedure
  # tolLasso : tolerance for stopping criterion in Lasso minimization
  # PostLasso : if TRUE computes the PostLasso instead of Lasso
  # trace : if TRUE print convergence info
  
  # The first column of X must be the constant

  
  ### Load user-defined functions
  source("functions/LassoFISTA.R")
  source("functions/LeastSq.R")
  source("functions/LeastSqgrad.R")
  library('lbfgs')
  
  ### Setting
  d <- as.matrix(d)
  y <- as.matrix(y)
  X <- as.matrix(X)
  
  n <- nrow(X)
  p <- ncol(X)
  
  ### Type of method
  if(method == "WLSLasso"){ W=(1-d)*exp(X%*%beta) }
  if(method == "LinearOutcome"){ W=(1-d)*(sum(d)/sum(1-d)) }
  W <- as.vector(W)
  
  ### First step: Lasso
  
  # Overall penalty level
  g <- .1/log(max(p,n))
  lambda <- c*qnorm(1-.5*g/p)/sqrt(n)
  
  # Penalty loadings: get a preliminary estimate
  m_y <- as.vector( t(W)%*%y/sum(W) )
  Psi <- diag(as.vector(sqrt( t(W*(y-m_y)^2) %*% (diag(sqrt(W))%*%X)^2 / n )))
  
  # Estimation parameters
  v <- .01 # Stopping rule
  k <- 0
#  mu <- c(m_y,rep(0,p-1))
  mu <- rep(0,p)
  
    
    # Lasso estimate
    repeat{
      k <- k+1
      
      # Compute Lasso estimate
      # LassoEstim <- LassoFISTA(betaInit=Psi%*%mu,y=y,X%*%solve(Psi),W=W,
      #                        lambda,tol=tolLasso,maxIter=maxIterLasso,trace=F)
		  # mu <- solve(Psi) %*% LassoEstim$beta
		  
      LassoEstim <- lbfgs(LeastSq, LeastSqgrad, Psi%*%mu, y=y, X=X%*%solve(Psi), W=W,
		                  orthantwise_c=lambda,orthantwise_start=1,
		                  invisible=1)
		  mu <- solve(Psi) %*% LassoEstim$par
      

      
      # Update penalty loadings
      PrePsi <- Psi
      Psi <- diag(as.vector(sqrt( t((W*y-(diag(W)%*%X)%*%mu)^2) %*% (diag(W)%*%X)^2 / n )))
      
      # Trace showing
      if(trace & k%%5==0){
        	print(paste("Max. pen. loading diff at Lasso Iteration nb.",k,":",max(abs(diag(Psi-PrePsi))))) 
      }
      
      # Stopping rules
#      if(k > maxIterPen || max(abs(diag(Psi-PrePsi))) < v || LassoEstim$convergenceFISTA==-555){
#        cvg <- LassoEstim$convergenceFISTA
#        break
 #     } 

if(k > maxIterPen || max(abs(diag(Psi-PrePsi))) < v){
  cvg <- LassoEstim$convergenceFISTA
  break
} 

    }
  
    cvg =0
    if(k > maxIterPen){
      cvg=-999
      print("Penalty estimation did not converge.")
    }

    muLasso <- mu 
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
  return(list(method=method,
              lambda=lambda,
              SHat=SHat,
              muPL=muPL,
              muLasso=c(muLasso),
              nbIter=k,
              convergence=cvg
              ))
  
  
}