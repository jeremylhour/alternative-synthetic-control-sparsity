#' Function to compute theta, the Immunized ATT
#' 
#' Return t-stat based on confidence interval
#' 
#' Last edited: 18 avril 2016.
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
#' @return tstat T Statistics for testing null 'theta=0'.
#'
#' @author Jeremy Lhour

ImmunizedATTCI <- function(theta,y,d,X,beta,mu=rep(0,ncol(X)),Immunity=T){
  
  d <- as.matrix(d)
  y <- as.matrix(y)
  X <- as.matrix(X)
  
  eps <- y
  if(Immunity) eps <- y - X%*%mu 
  
  pi <- mean(d)
  theta_hat <- mean((d - (1-d)*exp(X%*%beta)) * eps)  / pi
  
  # Compute standard error
  psi <- (d - (1-d)*exp(X%*%beta)) * eps - d*theta
  VAR <- mean(psi^2)/pi^2
  standard_dev <- sqrt(VAR)/sqrt(nrow(X))
  
  # Compute t-statistics for H_0: 'theta=0'
  psi_0 <- (d - (1-d)*exp(X%*%beta)) * eps
  VAR_0 <- mean(psi_0^2)/pi^2
  tstat <- sqrt(nrow(X)) * theta / sqrt(VAR_0)
  
  return(list(tstat=(theta_hat-theta)/standard_dev))
}