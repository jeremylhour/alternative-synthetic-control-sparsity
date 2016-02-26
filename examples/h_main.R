### Point sur h.
### Jeremy L Hour
### 26 fevrier 2016

rm(list=ls())

### Loi Normale
G_0 <- function(x){
  pnorm(x)
} 

G_1 <- function(x){
  dnorm(x)
}

G_2 <- function(x){
  -x*exp(-.5*x^2)/sqrt(2*pi)
}

### Loi de Cauchy
G_0 <- function(x){
  pcauchy(x)
} 

G_1 <- function(x){
  dcauchy(x)
}

G_2 <- function(x){
  -(2/pi)*x/(x^2+1)^2
}


### Define h
h_0 <- function(x){
  G_0(x)/(1-G_0(x))
}

h_1 <- function(x){
  G_1(x) /(1-G_0(x))^2
}

h_2 <- function(x){
  G_2(x)/(1-G_0(x))^2 + 2*G_1(x)^2/(1-G_0(x))^3
}

### 
x_val <- seq(-3,3,by=.1)

plot(x_val,h_0(x_val), type="line")
plot(x_val,h_1(x_val), type="line")
plot(x_val,h_2(x_val), type="line")