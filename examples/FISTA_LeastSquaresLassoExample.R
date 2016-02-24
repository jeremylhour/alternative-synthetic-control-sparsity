### Lasso minimization : Nesterov example
### Jeremy L Hour
### 14 janvier 2016

### Set working directory
setwd("R:/Simulations/R_Code")

rm(list=ls())
set.seed(9081993)

### Load packages
library("MASS")
library("lars")
library("lbfgs")

### Load user-defined functions
source("functions/LassoFISTA.R")
source("functions/DataSim.R") 
source("functions/LeastSq.R")
source("functions/LeastSqgrad.R")


### Select penalty level by cross-validation

data <- DataSim(n=500,p=400)
X <- data$X
y <- data$y
n <- nrow(X)
p <- ncol(X)

# Overall penalty level
LassoPen <- function(c=1.1,n,p){
  g <- .1/log(max(p,n))
  lambda <- c*qnorm(1-.5*g/p)/sqrt(n)
  return(lambda)
}

# 10-fold cross-validation
c_set <- seq(0,1.1,by=0.05)
fold <- sample(rep(1:10,each=n/10))
beta0 <- rep(0,p)

cv_error <- matrix(NA,nrow=10,ncol=length(c_set))
j <- 0

for(c in c_set){
  j <- j+1
  print(paste("penalty level:",c))
  for(v in 1:10){
    # Compute model
    X_s2 <- X[which(fold!=v),] 
    y_s2 <- y[which(fold!=v)] 
    fit <- LassoFISTA(betaInit=beta0,y_s2,X_s2,
               lambda=LassoPen(c=c,n=nrow(X_s2),p=ncol(X_s2)),trace=F,maxIter=10e6)
    
    # Estimate on the left-out sample
    X_s1 <- X[which(fold==v),] 
    y_s1 <- y[which(fold==v)] 
    cv_error[v,j] <- LeastSq(fit$beta,y_s1,X_s1,W=rep(1,nrow(X_s1)))
  }
}

# Best value of C by cross-validation
c_set[which(apply(cv_error,2,mean)== min(apply(cv_error,2,mean)))]



### MC XP
R <- 100
Results <- matrix(ncol=3, nrow=R)

for(r in 1:R){
  ### 1. Generate data
  print(paste("MC Iteration nb.",r))
  
  data <- DataSim(n=500,p=400)
  X <- data$X
  y <- data$y
  d <- data$d
  n <- nrow(X)
  p <- ncol(X)

beta0 <- rep(0,p)

# Overall penalty level
c <- 1.1
g <- .1/log(max(p,n))
lambda <- c*qnorm(1-.5*g/p)/sqrt(n)


# 2. FISTA
FISTAMethod <- LassoFISTA(betaInit=beta0,y,X,lambda=lambda,trace=F,maxIter=10e6)

# 3. LBFGS
LBFGS <- lbfgs(LeastSq, LeastSqgrad, beta0, y=y, X=X, W=rep(1,nrow(X)),
               invisible=1, orthantwise_c=lambda,
               linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",
               orthantwise_start = 1)
# 4. LARS
LARS <- lars(as.matrix(X),y, type="lasso", normalize=F, intercept=F,use.Gram=T)
betaLARS <- coef(LARS,s=sum(FISTAMethod$beta!=0))

Results[r,] <- c(LeastSq(FISTAMethod$beta,y,X,W=rep(1,nrow(X))) + lambda*sum(abs(FISTAMethod$beta)),
                 LeastSq(LBFGS$par,y,X,W=rep(1,nrow(X))) + lambda*sum(abs(LBFGS$par)),
                 LeastSq(betaLARS,y,X,W=rep(1,nrow(X))) + lambda*sum(abs(betaLARS)))
}