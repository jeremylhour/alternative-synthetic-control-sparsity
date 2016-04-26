### Trying stuff:
### Gradient of the estimating moment
psi <- function(y,d,X,beta,mu=rep(0,ncol(X)),theta){
  d <- as.matrix(d)
  y <- as.matrix(y)
  X <- as.matrix(X)
  
  w <- (d-(1-d)*exp(X%*%beta))*(y - X%*%mu) - d*theta
  
  return(mean(w))
}

dpsi_beta <- function(y,d,X,beta,mu=rep(0,ncol(X))){
  d <- as.matrix(d)
  y <- as.matrix(y)
  X <- as.matrix(X)
  
  w <- (d-1)*exp(X%*%beta)*(y - X%*%mu)
  return(t(w)%*%X/nrow(X))
}

dpsi_mu <- function(y,d,X,beta,mu=rep(0,ncol(X))){
  d <- as.matrix(d)
  y <- as.matrix(y)
  X <- as.matrix(X)
  
  w <- -(d - (1-d)*exp(X%*%beta))
  return(t(w)%*%X/nrow(X))
}

# BEAST
theta_BEAST <- ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$theta
psi(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso,theta_BEAST)

# Taylor approx
mean(d)*sqrt(nrow(X))*theta_BEAST # lhs
sqrt(nrow(X)) * psi(y,d,X,c(0,data$b),c(0,data$g),0)

eta_hat <- c(CAL$betaLasso,ORT_WLS_L$muLasso)
eta0 <- c(0,data$b,0,data$g)
sqrt(nrow(X)) * t(eta_hat-eta0) %*% c(dpsi_beta(y,d,X,c(0,data$b),c(0,data$g)), dpsi_mu(y,d,X,c(0,data$b),c(0,data$g)))

sum(abs(dpsi_beta(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso)))
sum(abs(dpsi_beta(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL)))

sum(abs(dpsi_beta(y,d,X,CAL$betaLasso,rep(0,ncol(X)))))
sum(abs(dpsi_beta(y,d,X,CAL$betaPL,rep(0,ncol(X)))))

BCHd <- coef(lm(d ~ X))

# Farrell
theta_FARRELL <- ImmunizedATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso, Immunity=T)$theta
psi(y,d,X,LOGIT$betaLasso, FARRELL$muLasso,theta_FARRELL)

sum(abs(dpsi_beta(y,d,X,LOGIT$betaLasso, FARRELL$muLasso)))
sum(abs(dpsi_beta(y,d,X,LOGIT$betaPL, FARRELL$muPL)))

sum(abs(dpsi_beta(y,d,X,LOGIT$betaLasso, rep(0,ncol(X)))))
sum(abs(dpsi_beta(y,d,X,LOGIT$betaPL, rep(0,ncol(X)))))