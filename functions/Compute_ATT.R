#' Function to compute the average treatment effect on the treated (ATT)
#' works for immunized, naive and low-dimensional.
#' compute immunized only if 'mu' is specified.
#' 
#' Replaces "ImmunizedATT", "LowDim_ATT" and "Naive_ATT"
#' 
#' Created: 13/05/2020
#' 
#' @param y Outcome variable
#' @param d Treatment indicator
#' @param X Matrix of covariates
#' @param beta Estimated Calibration parameter (step 1), can be sparse or not
#' @param mu Estimated Orthogonality parameter (step 2), optionnal. If specified return the immunized estimator
#' 
#' @return theta ATT estimate
#' @return sigma Asymptotic standard error
#' @return tstat T-statistics for testing null 'theta = 0'
#'
#' @author Jeremy L'Hour

Compute_ATT <- function(y,d,X,beta,mu){
  d <- as.matrix(d); y <- as.matrix(y); X <- as.matrix(X)
  pi <- mean(d)
  
  if(missing(mu)){
    theta = mean((d - (1-d)*exp(X%*%beta)) * y)  / pi # not immunized
    S_hat = which(beta != 0) # if low-dimensional, all should be different from zero.
    untreated_reg = lm(y ~ X[,S_hat] - 1, weights=(1-d)*exp(X%*%beta))
    mu = coef(untreated_reg)
    mu[is.na(mu)] = 0
    eps = y - X[,S_hat]%*%as.matrix(mu)
  } else {
    eps = y - X%*%mu
    theta = mean((d - (1-d)*exp(X%*%beta)) * eps)  / pi # immunized
  }
  
  # Compute standard error
  psi = (d - (1-d)*exp(X%*%beta)) * eps - d*theta
  VAR = mean(psi^2)/pi^2
  standard_dev = sqrt(VAR)/sqrt(nrow(X))
  
  # Compute t-statistics for H_0: 'theta=0'
  psi_0 = (d - (1-d)*exp(X%*%beta)) * eps
  VAR_0 = mean(psi_0^2)/pi^2
  tstat = sqrt(nrow(X)) * theta / sqrt(VAR_0)
  
  return(list(theta=theta,
              sigma=standard_dev,
              tstat=tstat))
}