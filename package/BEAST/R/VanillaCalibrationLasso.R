VanillaCalibrationLasso <- function(d,X,c=1.1,
                            PostLasso=F,trace=F,maxIter=1e3){
  ### Function to compute beta_0, the covariate balancing weights
  ### Only use a single penalty parameter
  ### Jeremy L Hour
  ### 11 janvier 2016
  
  ### INPUTS:
  # d : Treatment indicator
  # X : Covariates
  # c : constant for the overall penalty level
  # trace : if TRUE print convergence info
  # PostLasso : if TRUE computes the PostLasso instead of Lasso
  # maxIter : maximal number of iteration for Lasso
  
  ### Load necessary packages
  library("lbfgs")
  
  ### Setting
  d <- as.matrix(d)
  X <- as.matrix(X)
  
  n <- nrow(X)
  p <- ncol(X)
  
  ### First step: Lasso
  
  # Overall penalty level
  g <- .1/log(max(p,n))
  lambda <- c*qnorm(1-.5*g/p)/sqrt(n)
  
  
  beta <- c(log(sum(d)/(sum(1-d))),rep(0,p-1))
  # Compute Lasso estimate, do not penalize constant!
  LassoEstim <- lbfgs(gamma, gammagrad, beta, d=d, X=X,
                        orthantwise_c=lambda,orthantwise_start=1,
                        invisible=1)
    
  # Get coefficient
  betaLasso <- LassoEstim$par

  # Get estimated active set
  SHat <- union(1,which(betaLasso != 0)) #Always put a constant
  
  
  ### Second step: Post-Lasso
  betaPL <- rep(0,p)
  if(PostLasso==T){
    PL <- lbfgs(gamma, gammagrad, beta[SHat], d=d, X=X[,SHat],
                orthantwise_c=0, invisible=1)
    betaPL[SHat] <- PL$par
  }
  
  
  
return(list(SHat=SHat,
            betaPL=betaPL,
            betaLasso=c(betaLasso),
            lambda=lambda
            ))
  
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