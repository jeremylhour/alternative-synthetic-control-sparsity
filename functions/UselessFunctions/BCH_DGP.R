### BCH "Inference after model selection", ReStud 2013
### DGP for Monte Carlo experiment
### 4 fevrier 2016
### J L'Hour


BCH_DGP <- function(n=100,p=200,cy=.824,cd=1, Design=1){
  
    a <- 0.5 # Treatment effect
  
    mu <- rep(0,p)
    
    for(j in 1:p){
      mu[j] <- 1 / j^2
    }
    
    
  ### Treatment variable coefficient
  b <- cy * mu
  ### Outcome equation coefficients
  gamma <- cd * mu
    
  ### Correlation coefficients between x's
  rho <- .5 # correlation coefficient
  Sigma <- matrix(0,nrow=p, ncol=p)
  
  for(k in 1:p){
    for(j in 1:p){
      Sigma[k,j] <- rho^abs(k-j)
    }
  }
  
  # Simulate covariates
  X <- mvrnorm(n = n, mu=rep(0,p), Sigma)
  
  
  if(Design==1){
    # Simulate treatment
    d <- as.numeric(X%*%gamma + rnorm(n) > 0)

    # Simulate outcome
    y <- a*d + X%*%b + rnorm(n)
  }
  
  
  if(Design==2){
    ### Treatment error variance (Design 2)
    sd <- sqrt((1 + X%*%mu)^2 / mean((1 + X%*%mu)^2))
    
    # Simulate treatment
    d <- as.numeric(X%*%gamma +  sd*rnorm(n) > 0)
    
    ### Outcome error variance (Design 2)
    sy <- sqrt((1 + a*d + X%*%mu)^2 / mean((1 + a*d + X%*%mu)^2))
    
    # Simulate outcome
    y <- a*d + X%*%b + sy*rnorm(n)
  }
  
  
 
  
  # Add the intercept
  X <- cbind(rep(1,n),X)
  
  return(list(X=X,y=y,d=d, b=b,g=gamma))
}