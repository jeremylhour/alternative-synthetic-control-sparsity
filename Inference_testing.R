### High-dimension Calibration for Treatment Effect
### Inference and testing
### 18 avril 2015
### J L'Hour

### Set working directory
setwd("R:/Simulations/BEAST") # R pour Jeremy, Z pour Marianne

rm(list=ls())
set.seed(30031987)


### 0. Settings

### Load packages
library("ggplot2")
library("gridExtra")

### Load user-defined functions
source("functions/DataSim.R") 
source("functions/DataSimCauchy.R") 
source("functions/AwkwardDataSim.R") 
source("functions/CalibrationLasso.R")
source("functions/OrthogonalityReg.R")
source("functions/LogitLasso.R")
source("functions/BCHDoubleSelec.R")
source("functions/VanillaCalibrationLasso.R")
source("functions/VanillaOrthogonalityReg.R")
source("functions/ImmunizedATT.R")


### MC XP FUNCTION
R <- 500

MonteCarloSimu <- function(Ry=.8,Rd=.2, n=250,p=300){
  Results <- matrix(ncol=1, nrow=R)
  t_start <- Sys.time()
  pb <- txtProgressBar(style = 3)
  
  for(r in 1:R){
    ### 1. Generate data
    data <- DataSim(n=n,p=p,Ry=Ry,Rd=Rd,TreatHeter=F)
    X <- data$X
    y <- data$y
    d <- data$d
    
    ### 2. Calibration part
    CAL <- CalibrationLasso(d,X,c=.5,maxIterPen=1,PostLasso=F,trace=F,maxIter=1e6)
    
    ### 3. Computes the orthogonality parameter, using method WLS Lasso
    ORT_WLS_L <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="WLSLasso",
                                  c=2, nopenset=c(1), RescaleY=F,
                                  maxIterPen=1,maxIterLasso=1e4,tolLasso=1e-6,PostLasso=F,trace=F)
    
    ### 3. Third step: ATT estimation
    Results[r,] <- ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$tstat
    
    setTxtProgressBar(pb, r/R)
    
  }
  
  close(pb)
  print(Sys.time()-t_start)
  return(Results)
}

### RUN SIMULATIONS
i <- 0
j <- 0
z <- matrix(,ncol=10,nrow=10)

for(Ry in seq(0,9,by=.1)){
  i <- i+1
  print(paste("Ry:",i))
  for(Rd in seq(0,9,by=.1)){
    j <- j+1
    print(paste("Rd:",j))
    MCX <- MonteCarloSimu(Ry=Ry,Rd=Rd, n=250,p=300)
    z[i,j] <- sum(ifelse(abs(MCX) > qnorm(.975),1,0))/R
  }
}


# H_0 should be rejected 5 pcts of the time
sum(ifelse(abs(Results) > qnorm(.975),1,0))/R

# Draw the charts
id <- c(mapply(function(x) rep(x,R),1))
val <- c(Results)
data_res <- data.frame(val = val, model = id)

M <- max(abs(quantile(Results,.01)),abs(quantile(Results,.99)))
lb <- -1.1*M
ub <- 1.1*M


### Function for plot
get.plot <- function(data,modelS,title="A Title",sdBCH){
  plot_res <- ggplot(subset(data, (model==modelS)), aes(x=val)) + 
    geom_histogram(binwidth = .2, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
    scale_x_continuous(limits=c(lb,ub), name="Treatment effect") +
    ggtitle(title) + 
    stat_function(fun = dnorm, args=list(mean=0, sd=sdBCH), colour="darkorchid3", size=1) +
    theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")
  
  return(plot_res)
}


get.plot(data_res,1,"T-Stat Distribution", 1)



### Function

g <- function(x, theta){
  return( (theta/(1+theta^2)) * exp(-0.5*(x-theta)^2) )
}

x=seq(-10,10,by=.5)
theta=seq(-10,10,by=.5)
z=outer(x,theta,g)

op <- par(bg = "white")
pdf("Confidence_covering_shape.pdf")
persp(x,theta,z,
      theta=30,
      phi=20,
      expand=0.5,
      col = "steelblue",
)
dev.off()