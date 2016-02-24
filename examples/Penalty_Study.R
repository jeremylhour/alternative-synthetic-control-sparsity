lambda <- function(n,p,c=1.1){
  g <- .1/log(max(p,n))
  lambda <- c*qnorm(1-.5*g/p)/sqrt(n)
  return(lambda)
}

n <- seq(100,10000,by=100)
p <- seq(10,500,by=50)
my.colors<- rainbow(length(p))

plot(n,lambda(n,p[1]), type="line",lwd=2, ylim=c(0,.4),col=my.colors[1])

for(j in 2:length(p)){
  lines(n,lambda(n,p[j]),lwd=2, col=my.colors[j])
}

### rationaliser le c de l'application tabac 
lambda_applitabac <- 6.4113e-04
n <- 39
p <- 13

k <- 0
c <- 1.01006e-05
repeat{
  k <- k+1
  if(lambda(n,p,c)-lambda_applitabac > 0){
    c <- .99999*c
  } else {
    c <- 1.00001*c
  }
  if(k > 1000) break
}
