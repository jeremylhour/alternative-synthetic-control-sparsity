#' Function to compute theta, ATT in low-dimension
#' 
#' Last edited:12 mai 2020.
#' 
#' @param y Outcome variable.
#' @param d Treatment indicator.
#' @param X Matrix of covariates.
#' @param beta Estimated Calibration parameter (Step 1).
#' 
#' @return theta Immunized ATT estimate.
#' @return sigma Asymptotic standerd error. 
#' @return tstat T Statistics for testing null 'theta=0'.
#'
#' @author Jeremy L'Hour

LowDim_ATT <- function(y,d,X,beta){
  d = as.matrix(d); y = as.matrix(y); X = as.matrix(X)
  
  OrthoReg = lm(y ~ X - 1, weights=(1-d)*exp(X%*%beta))
  mu = coef(OrthoReg)
  
  pi = mean(d)
  theta = mean((d - (1-d)*exp(X%*%beta)) * y)  / pi
  
  # Compute standard error
  psi = (d - (1-d)*exp(X%*%beta)) * (y - X%*%mu) - d*theta
  VAR = mean(psi^2)/pi^2
  standard_dev = sqrt(VAR)/sqrt(nrow(X))
  
  # Compute t-statistics for H_0: 'theta=0'
  psi_0 = (d - (1-d)*exp(X%*%beta)) * (y - X%*%mu)
  VAR_0 = mean(psi_0^2)/pi^2
  tstat = sqrt(nrow(X)) * theta / sqrt(VAR_0)
  
  return(list(theta=theta,
              sigma=standard_dev,
              tstat=tstat))
}