#' New DGP for Monte-Carlo Experiment
#' 
#' Created: May 7, 2020
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

New_DataSim <- function(n=2000, p=50, Ry=.5, Rd=.2, Intercept=T, rho=.5){
  # Variance matrix for X
  Sigma = matrix(0,nrow=p, ncol=p)
  for(k in 1:p){
    for(j in 1:p){
      Sigma[k,j] = rho^abs(k-j)
    }
  }
  
  # Treatment selection coefficient
  gamma = rep(0,p)
  for(j in 1:10){
    #gamma[j] = 1*(-1)^(j) / j^2
    gamma[j] = 1 / j^2
  }
  
  # Outcome equation coefficients
  mu = gamma
  for(j in (p-9):p){
    mu[j] = (-1)^(j+1) / (p-j+1)^2
  }
  
  # Adjustment to match R.squared in latent equation
  c_gamma = c(sqrt((1/t(gamma)%*%Sigma%*%gamma)*(Rd/(1-Rd))))
  gamma = c_gamma*gamma
  
  # Adjustment to match R.squared of Y_0
  c_mu = c((Ry/(2*(1-Ry)*(t(mu)%*%Sigma%*%mu)^4))^(1/8))
  
  # Simulating the DGP
  X = mvrnorm(n=n, mu=rep(0,p), Sigma)

  
  d = as.numeric(runif(n) < sigmoid(X%*%gamma))
  
  y = (X%*%(c_mu*mu))^2 + d*(X%*%gamma) + rnorm(n)
  
  # Compute the ATT
  Z = rnorm(5000000, sd=t(gamma)%*%Sigma%*%gamma)
  ATT = mean(Z*sigmoid(Z)/c_gamma)/mean(sigmoid(Z))
  
  if(Intercept) X <- cbind(rep(1,n),X)
  
  return(list(X=X,
              y=y,
              d=d,
              mu=mu,
              gamma=gamma,
              ATT=ATT))
}

sigmoid = function(x) 1/(1+exp(-x))