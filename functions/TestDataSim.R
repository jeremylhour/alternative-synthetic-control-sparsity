#' New DGP for Monte-Carlo Experiment
#' 
#' Created: May 19, 2021
#' 
#' @param n Sample size
#' @param p Number of Covariates, must be larger than 10
#' @param Ry Desired R-squared for outcome equation
#' @param Rd Desired R-squared for treatment equation
#' @param Intercept If TRUE returns an intercept also
#' @param rho Correlation coefficients between two adjacent covariates.
#' 
#' @return Returns y, d and X necessary data to run Monte Carlo simulations. Also returns parameters value for a possible check.
#' 
#' @author Jeremy L'Hour

TestDataSim <- function(n=2000, p=50, Ry=.5, Rd=.2, Intercept=T, rho=.5){
  # Variance matrix for X
  Sigma = matrix(0,nrow=p, ncol=p)
  for(k in 1:p){
    for(j in 1:p){
      Sigma[k,j] = rho^abs(k-j)
    }
  }
  
  # Treatment selection coefficient
  gamma_0 = rep(0,p)
  for(j in 1:10){
    gamma_0[j] = 1*(-1)^(j) / j^2
  }
  
  # Outcome equation coefficients
  mu =  gamma_0
  for(j in (p-9):p){
    mu[j] = (-1)^(j+1) / (p-j+1)^2
  }
  
  # Adjustment to match R.squared in latent equation
  c_gamma = c(sqrt((1/t(gamma_0)%*%Sigma%*%gamma_0)*(Rd/(1-Rd))))
  gamma_0 = c_gamma*gamma_0
  
  # Adjustment to match R.squared of Y_0
  c_mu = sqrt(log(Ry/(1-Ry))/(2*t(mu)%*%Sigma%*%mu))
  mu = c(c_mu)*mu
  
  # Adjustment of heterogeneity
  eta = sqrt( (exp(t(mu)%*%Sigma%*%mu) -1)*exp(t(mu)%*%Sigma%*%mu)/(4*t(gamma_0)%*%Sigma%*%gamma_0) )
  
  # Simulating the DGP
  X = mvrnorm(n=n, mu=rep(0,p), Sigma)
  
  d = as.numeric(runif(n) < link_func(X%*%gamma_0))
  
  y = exp(X%*%mu) + c(eta)*d*(X%*%gamma_0) + rnorm(n)
  
  # Compute the ATT
  Z = rnorm(100000, sd=c(eta)*sqrt(t(gamma_0)%*%Sigma%*%gamma_0))
  ATT = mean(c(eta)*Z*link_func(Z))/mean(link_func(Z))
  
  if(Intercept) X <- cbind(rep(1,n),X)
  
  return(list(X=X,
              y=y,
              d=d,
              mu=mu,
              gamma_0=gamma_0,
              ATT=ATT))
}

link_func = function(x) x^2/(1+x^2)