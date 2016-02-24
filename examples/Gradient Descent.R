### Lasso avec Nesterov

### See http://www.seas.ucla.edu/~vandenbe/ee236c.html

### With our function

# Overall penalty level
g <- .1/log(max(p,n))
lambda <- c*qnorm(1-.5*g/p)/sqrt(n)

# Penalty loadings: get a preliminary estimator
beta <- 0
prelim <- lbfgs(gamma, gammagrad, beta, Z=cbind(d,X[,1]),
                invisible=1, orthantwise_c=0,
                linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",
                orthantwise_start = 1)
beta <- c(prelim$par,rep(0,p-1))
Psi <- diag(as.vector(sqrt( t((d - exp(X%*%beta)*(1-d))^2) %*% X^2 / n )))

# Compute Lasso estimate
LassoEstim <- lbfgs(gamma, gammagrad, beta, Z=cbind(d,X %*% solve(Psi)),
                    invisible=1, orthantwise_c=lambda,
                    linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",
                    orthantwise_start = 1)

# Get coefficient
beta <- solve(Psi) %*% LassoEstim$par


### 2. With Nesterov

t <- .001
K <- 5000
beta <- 0
prelim <- lbfgs(gamma, gammagrad, beta, Z=cbind(d,X[,1]),
                invisible=1, orthantwise_c=0,
                linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",
                orthantwise_start = 1)
beta <- c(prelim$par,rep(0,p-1))
Psi <- diag(as.vector(sqrt( t((d - exp(X%*%beta)*(1-d))^2) %*% X^2 / n )))
v <- beta

k <- 0
tol <- 1e-6

repeat{
  k <- k+1
  print(k)
  print(paste("Difference:",abs(gamma(beta,cbind(d,X %*% solve(Psi)))-gamma(betaO,cbind(d,X %*% solve(Psi))))))
  
  theta <- 2/(k+1)
  betaO <- beta
  z <- (1-theta)*beta + theta * v
  v <- prox(v - (t/theta) * gammagrad(z,cbind(d,X %*% solve(Psi))), lambda*t/theta)
  beta <- (1-theta)*beta + theta * v
  
  if(abs(gamma(beta,cbind(d,X %*% solve(Psi)))-gamma(betaO,cbind(d,X %*% solve(Psi)))) < tol || k > K) break 
}

# Get coefficient
beta <- solve(Psi) %*% LassoEstim$par







### Define appropriate functions
MC <- function(beta,y,X){
  n <- length(y)
  return(t(y-X%*%beta) %*% (y-X%*%beta) / n)
}

MCgrad <- function(beta,y,X){
  n <- length(y)
  return( -2 * t(X) %*% (y-X%*%beta)/ n)
}

prox <- function(x,lambda){
  (abs(x)-lambda)*( abs(x) -lambda > 0) * sign(x)
}


## Set starting point
lambda <- .5 # Lasso penalty
t <- .5*n/max(eigen(t(X)%*%X)$values)
K <- 5000
beta <- rep(0,p)
v <- beta

k <- 0

repeat{
  k <- k+1
  print(k)
  
  theta <- 2/(k+1)
  z <- (1-theta)*beta + theta * v
  v <- prox(v - (t/theta) * MCgrad(z,y,X), lambda*t/theta)
  beta <- (1-theta)*beta + theta * v
  
  if(k > K) break 
}

LassoEstim <- lbfgs(MC, MCgrad, beta, y=y, X=X,
                    invisible=1, orthantwise_c=lambda,
                    linesearch_algorithm="LBFGS_LINESEARCH_BACKTRACKING",
                    orthantwise_start = 1)