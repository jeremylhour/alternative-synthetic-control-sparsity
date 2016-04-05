
LassoFISTADebug <- function(betaInit=rep(0,ncol(X)),y,X,W=rep(1,nrow(X)),
                       nopen=NULL,lambda, psi=rep(1,ncol(X)),
                       tol=1e-8,maxIter=1000,trace=F){
  # Observation weighting
  W <- as.vector(W)
  y <- sqrt(W)*y
  X <- diag(sqrt(W)) %*% as.matrix(X)
  
  # Penalty loadings
  if(sum(psi<0)>1) stop("Can only use non-negative penalty loadings!")
  psi <- as.vector(psi)
  if(length(psi)!=ncol(X)) stop("Size of penalty loadings vector must be of ncol(X).")
  
  ### Set Algo. Values
  eta <- 1/max(2*eigen(t(X)%*%X)$values/nrow(X))
  theta <- 1
  thetaO <- theta
  beta <- betaInit
  v <- beta
  cv <- 0
  
  f_path <- NULL
  
  k <- 0
  repeat{
    k <- k+1
    
    thetaO <- theta
    theta <- (1+sqrt(1+4*thetaO^2))/2
    delta <- (1-thetaO)/theta
    
    betaO <- beta
    beta <- prox(v - eta*LeastSqgrad(v,y,X), lambda*eta, psi, nopen)
    
    v <- (1-delta)*beta + delta*betaO
    
    f_path <- c(f_path,LassoObj(beta,y,X,lambda,psi,nopen))
    
    # Show objective function value
    if(trace==T & k%%100 == 0){ print(paste("Objective Func. Value at iteration",k,":",LassoObj(beta,y,X,lambda,psi,nopen))) }
    
    # Break if diverges
    if(is.na(LassoObj(beta,y,X,lambda,psi,nopen) - LassoObj(betaO,y,X,lambda,psi,nopen))){
      cv <- -555
      print("LassoFISTA did not converge")
      break
    } else if(abs(LassoObj(beta,y,X,lambda,psi,nopen) - LassoObj(betaO,y,X,lambda,psi,nopen)) < tol || k > maxIter) break
    
  }
  
  if(k > maxIter){
    print("Reach max. number of iterations reach in Lasso minimization.")
    cv <- -555
  } 
  
  return(list(beta=beta,
              value=LassoObj(beta,y,X,lambda,psi,nopen),
              loss=LeastSq(beta,y,X),
              l1norm=abs(beta),
              nbIter=k,
              convergenceFISTA=cv,
              f_path=f_path))
}


#################################
#################################
### Define auxiliary functions###
#################################
#################################

prox <- function(x,lambda,psi,nopen){
  y <- (abs(x)-lambda)*(abs(x)-lambda*psi > 0) * sign(x)
  y[nopen] <- x[nopen] # Do not penalize these variables
  return(y)
}

LeastSq <- function(mu,y,X){
  X <- as.matrix(X)
  return(mean((y - X%*%mu)^2))
}

LeastSqgrad <- function(mu,y,X){
  X <- as.matrix(X)
  df <- as.vector(-2*(t(y - X%*%mu)%*%X) / nrow(X))
  return(df)
}

LassoObj <- function(beta,y,X,lambda,psi,nopen){
  if(length(nopen)>0){
    psi[nopen]=0
    f <- LeastSq(beta,y,X) + lambda*sum(psi*abs(beta))
  } else {
    f <- LeastSq(beta,y,X) + lambda*sum(psi*abs(beta))
  }
  return(f)
}


# Compute Lasso estimate
LassoEstim <- LassoFISTADebug(betaInit=mud,d,X,W=rep(1,n),
                         nopen=nopenset,lambda=lambda,psi=Psi,
                         tol=.001,maxIter=6000,trace=T)

fp <- LassoEstim$f_path

fp <- (fp-min(fp))/min(fp)

plot(fp[0:35], type="line")