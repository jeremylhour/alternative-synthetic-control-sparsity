#' Farrell (2015) Variance
#' 
#' Variance computation from Farrell (2015)
#' Last edited: 26 fevrier 2016.
#' 
#' @param y Outcome variable, not normalized.
#' @param X Matrix of covariates.  
#' @param d Treatment indicator
#' @param beta Logit coefficients
#' @param mu Outcome coefficient
#' 
#' @return sigma Variance for the ATT estimator
#' 
#' @author Jeremy Lhour


VarianceFarrell(y,X,d,beta,mu){
  source("functions/ImmunizedATT.R")
  
  pi <- mean(d)
  theta <- ImmunizedATT(y,d,X,beta,mu, Immunity=T)$theta
  eps <- y - X%*%mu
  
  va <- eps^2*(1-d)*(2/(1-pi)^2 + exp(2*X%*%beta)/pi - 2*exp(X%*%beta)/(pi*(1-pi))) + (eps-theta)^2*d/pi^2
  
  return(sqrt(mean(va)/nrow(X)))
  
}
