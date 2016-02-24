### Calibration cross-validation
### 22 fevrier 2016
### J L'Hour

### Set working directory
setwd("R:/Simulations/R_Code")

rm(list=ls())
set.seed(30031987)


### 0. Settings

### Load packages
library("ggplot2")
library("gridExtra")
library("MASS")
library("lbfgs")

### Load user-defined functions
source("functions/DataSim.R") 
source("functions/gamma.R")
source("functions/gammagrad.R")

### Data simulation
data <- DataSim(n=2000,p=100,Ry=.8,Rd=.2)
X <- scale(data$X)
X[,1] <- rep(1,nrow(X))
y <- data$y
d <- data$d

n <- nrow(X)
p <- ncol(X)

# Overall penalty level
LassoPen <- function(c=1.1,n,p){
  g <- .1/log(max(p,n))
  lambda <- c*qnorm(1-.5*g/p)/sqrt(n)
  return(lambda)
}

# 10-fold cross-validation
c_set <- seq(0,5,by=0.05)
fold <- sample(rep(1:10,each=n/10))
beta0 <- rep(0,p)

cv_error <- matrix(NA,nrow=10,ncol=length(c_set))
j <- 0
t_start <- Sys.time()
beta0 <- c(log(sum(d)/(sum(1-d))),rep(0,p-1))

for(c in c_set){
  j <- j+1
  print(paste("penalty level:",c))
  for(v in 1:10){
    # Compute model
    X_s2 <- X[which(fold!=v),] 
    d_s2 <- d[which(fold!=v)] 
    fit <-  lbfgs(gamma, gammagrad, beta0, d=d_s2, X=X_s2,
                               orthantwise_c=LassoPen(c=c,n=nrow(X_s2),p=ncol(X_s2)),
                              orthantwise_start=1,invisible=1)
    
    # Estimate on the left-out sample
    X_s1 <- X[which(fold==v),] 
    d_s1 <- d[which(fold==v)] 
    cv_error[v,j] <- gamma(fit$par,d_s1,X_s1)
  }
}
print(Sys.time()-t_start)


# Best value of C by cross-validation
mcv <- c_set[which(apply(cv_error,2,sum)== min(apply(cv_error,2,sum)))]

# Plot the result
plot(c_set,apply(cv_error,2,sum),
     col="steelblue", pch=20,
     ylim=c(5,7),
     xlab="Penalty value",
     ylab="Cross-validation error")
abline(v=mcv,
       lty=2,col="firebrick")


####################
####################
####################
## LALONDE DATASET #
####################
####################
####################

load("datasets/LalondeData_Unscaled.R")
X <- data$X
d <- data$d

n <- nrow(X)
p <- ncol(X)

# 10-fold cross-validation
c_set <- seq(0,4,by=0.05)
fold <- sample(rep(1:10,each=n/10))
beta0 <- rep(0,p)

cv_error <- matrix(NA,nrow=10,ncol=length(c_set))
j <- 0
t_start <- Sys.time()

beta0 <- lbfgs(gamma, gammagrad, beta0, d=d, X=X,
               orthantwise_c=.16,
               orthantwise_start=1,invisible=0)
beta0 <- beta0$par

for(c in c_set){
  j <- j+1
  print(paste("penalty level:",c))
  for(v in 1:10){
    # Compute model
    X_s2 <- X[which(fold!=v),] 
    d_s2 <- d[which(fold!=v)] 
    fit <-  lbfgs(gamma, gammagrad, beta0, d=d_s2, X=X_s2,
                  orthantwise_c=LassoPen(c=c,n=nrow(X_s2),p=ncol(X_s2)),
                  orthantwise_start=1,invisible=1)
    beta0 <- fit$par
    
    # Estimate on the left-out sample
    X_s1 <- X[which(fold==v),] 
    d_s1 <- d[which(fold==v)] 
    cv_error[v,j] <- gamma(fit$par,d_s1,X_s1)
  }
}
print(Sys.time()-t_start)


# Best value of C by cross-validation
mcv <- c_set[which(apply(cv_error,2,sum)== min(apply(cv_error,2,sum)))]

# Plot the result
plot(c_set,apply(cv_error,2,sum),
     col="steelblue", pch=20,
     ylim=c(0,3),
     xlab="Penalty value",
     ylab="Cross-validation error")
abline(v=mcv,
       lty=2,col="firebrick")

  
  
  