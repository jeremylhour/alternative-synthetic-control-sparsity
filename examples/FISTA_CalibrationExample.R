### Calibration minimization: examples
### Jeremy L Hour
### 14 janvier 2016

### Set working directory
setwd("R:/Simulations/BEAST")

rm(list=ls())
set.seed(9081993)

### Load packages
library("MASS")
library("lbfgs")

### Load user-defined functions
source("functions/CalibrationFISTA.R")
source("functions/CalibrationFISTA_Descent.R")
source("functions/DataSim.R") 
source("functions/gamma.R") 
source("functions/gammagrad.R") 


### MC XP
R <- 100
Results <- matrix(ncol=4, nrow=R)
Times <- matrix(ncol=3, nrow=R)

for(r in 1:R){
  ### 1. Generate data
  print(paste("MC Iteration nb.",r))
  
  data <- DataSim(n=5000,p=50,Scenario="AS")
  X <- data$X
  y <- data$y
  d <- data$d
  n <- nrow(X)
  p <- ncol(X)
  
  beta0 <- rep(0,p)
  
  # Overall penalty level
  c <- 1.01
  g <- .1/log(max(p,n))
  lambda <- c*qnorm(1-.5*g/p)/sqrt(n)
  
  # Optimization
  FISTA <- CalibrationFISTA(betaInit=beta0,d,X,lambda,
                              tol=1e-8,maxIter=10e6,maxLineSearch=10,trace=T)
  FISTA_D <- CalibrationFISTA_Descent(betaInit=beta0,d,X,lambda,
                            tol=1e-6,maxIter=5e2,maxLineSearch=10,trace=T)  
  
  
  BFGS <- lbfgs(gamma, gammagrad, beta0, d=d, X=X,orthantwise_c=lambda, invisible=1)

  # Record results  
  Results[r,] <- c(FISTA$value,
                   FISTA_D$value,
                   BFGS$value)
  Times[r,] <- c(FISTA$duration,
             FISTA_D$duration,
             FISTA_N$duration)
}

plot((FISTA$f_path-min(FISTA$f_path))/min(c(FISTA$f_path,FISTA_D$f_path)),type="line",col="blue",lty="dashed",lwd=2)
lines((FISTA_D$f_path-min(FISTA_D$f_path))/min(c(FISTA$f_path,FISTA_D$f_path)),col="red",lty="dotted",lwd=2)


Logit <- glm(d ~ X-1, family = "binomial")
beta0 <- Logit$coef

