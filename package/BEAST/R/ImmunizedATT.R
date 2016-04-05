#' Function to compute theta, the Immunized ATT
#' 
#' Last edited: 16 fevrier 2016.
#' 
#' @param y Outcome variable.
#' @param d Treatment indicator.
#' @param X Matrix of covariates.
#' @param beta Estimated Calibration parameter (Step 1).
#' @param mu Estimated Orthogonality parameter (Step 2).
#' @param Immunity If TRUE computes the immunized ATT, simple plug-in otherwise.
#' 
#' @return theta Immunized ATT estimate.
#' @return sigma Asymptotic standerd error. 
#'
#' @author Jeremy Lhour

ImmunizedATT <- function(y,d,X,beta,mu=rep(0,ncol(X)),Immunity=T){
  
  d <- as.matrix(d)
  y <- as.matrix(y)
  X <- as.matrix(X)
  
  eps <- y
  if(Immunity) eps <- y - X%*%mu 
  
  pi <- mean(d)
  theta <- mean((d - (1-d)*exp(X%*%beta)) * eps)  / pi
  
  psi <- (d - (1-d)*exp(X%*%beta)) * eps - d*theta
  VAR <- mean(psi^2)/pi^2
  standard_dev <- sqrt(VAR)/sqrt(nrow(X))
  
  return(list(theta=theta,
              sigma=standard_dev))
}