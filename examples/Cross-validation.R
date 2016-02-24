### Lasso cross-validation
### 22 fevrier 2016
### J L'Hour

### Set working directory
setwd("R:/Simulations/BEAST")

rm(list=ls())
set.seed(30031987)


### 0. Settings

### Load packages
library("ggplot2")
library("gridExtra")
library("MASS")

### Load user-defined functions
source("functions/DataSim.R") 
source("functions/LassoFISTA.R")

### Data simulation
data <- DataSim(n=200,p=100,Ry=.8,Rd=.2)
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

# Loss function
LeastSq <- function(mu,y,X){
  X <- as.matrix(X)
  return(mean((y - X%*%mu)^2))
}

### Comparison with LARS
fit <- LassoFISTA(betaInit=rep(0,p),y,X,
                  nopen=c(1),lambda=LassoPen(c=1.1,n=nrow(X),p=ncol(X)),
                  tol=1e-6,maxIter=1e4,trace=T)
library("lars")
lars_fit <- lars(X,y, type="lasso", normalize=F,intercept=F,trace=T)

betalars <- coef(lars_fit, s=n*LassoPen(c=1.1,n=nrow(X),p=ncol(X))/2, mode="lambda")


LeastSq(fit$beta,y,X) + LassoPen(c=1.1,n=nrow(X),p=ncol(X))*sum(abs(fit$beta[-1]))
LeastSq(betalars,y,X) + LassoPen(c=1.1,n=nrow(X),p=ncol(X))*sum(abs(betalars[-1]))



# 10-fold cross-validation
c_set <- seq(0,5,by=0.05)
fold <- sample(rep(1:10,each=n/10))
beta0 <- rep(0,p)

cv_error <- matrix(NA,nrow=10,ncol=length(c_set))
j <- 0
t_start <- Sys.time()

for(c in c_set){
  j <- j+1
  print(paste("penalty level:",c))
  for(v in 1:10){
    # Compute model
    X_s2 <- X[which(fold!=v),] 
    y_s2 <- y[which(fold!=v)] 
    fit <- LassoFISTA(betaInit=rep(0,p),y_s2,X_s2,
               nopen=c(1),lambda=LassoPen(c=c,n=nrow(X_s2),p=ncol(X_s2)),
               tol=1e-6,maxIter=1e4,trace=T)
    
    # Estimate on the left-out sample
    X_s1 <- X[which(fold==v),] 
    y_s1 <- y[which(fold==v)] 
    cv_error[v,j] <- LeastSq(fit$beta,y_s1,X_s1)
  }
}
print(Sys.time()-t_start)


# Best value of C by cross-validation
mcv <- c_set[which(apply(cv_error,2,mean)== min(apply(cv_error,2,mean)))]

# Plot the result
plot(c_set,apply(cv_error,2,mean),
     col="steelblue", pch=20,
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

load("datasets/LalondeData.R")
X <- scale(data$X)
X[,1] <- rep(1,nrow(X))
y <- (data$y - mean(data$y))/sd(data$y)
d <- data$d

n <- nrow(X)
p <- ncol(X)

# 10-fold cross-validation
c_set <- seq(0,5,by=0.05)
fold <- sample(rep(1:10,each=n/10))
beta0 <- rep(0,p)

cv_error <- matrix(NA,nrow=10,ncol=length(c_set))
j <- 0
t_start <- Sys.time()
beta0 <- rep(0,p)

for(c in c_set){
  j <- j+1
  print(paste("penalty level:",c))
  for(v in 1:10){
    # Compute model
    X_s2 <- X[which(fold!=v),] 
    y_s2 <- y[which(fold!=v)] 
    fit <- LassoFISTA(betaInit=beta0,y_s2,X_s2,
                      nopen=c(1),lambda=LassoPen(c=c,n=nrow(X_s2),p=ncol(X_s2)),
                      tol=1e-6,maxIter=1e4,trace=T)
    beta0 <- fit$beta
    
    # Estimate on the left-out sample
    X_s1 <- X[which(fold==v),] 
    y_s1 <- y[which(fold==v)] 
    cv_error[v,j] <- LeastSq(fit$beta,y_s1,X_s1)
  }
}
print(Sys.time()-t_start)


# Best value of C by cross-validation
mcv <- c_set[which(apply(cv_error,2,mean)== min(apply(cv_error,2,mean)))]

# Plot the result
plot(c_set,apply(cv_error,2,mean),
     col="steelblue", pch=20,
     xlab="Penalty value",
     ylim=c(0,.5),
     ylab="Cross-validation error")
abline(v=mcv,
       lty=2,col="firebrick")



for(r in 1:R){
  ### 1. Generate data
  print(paste("MC Iteration nb.",r))
  
  
  