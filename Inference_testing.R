### High-dimension Calibration for Treatment Effect
### Inference and testing
### 18 avril 2015
### J L'Hour

### Set working directory
setwd("R:/Simulations/BEAST") # R pour Jeremy, Z pour Marianne

rm(list=ls())
set.seed(12071990)


### 0. Settings

### Load packages
library("ggplot2")
library("gridExtra")
library("MASS")

### Load user-defined functions
source("functions/CalibrationLasso.R")
source("functions/OrthogonalityReg.R")
source("functions/LogitLasso.R")
source("functions/BCHDoubleSelec.R")
source("functions/VanillaCalibrationLasso.R")
source("functions/VanillaOrthogonalityReg.R")
source("functions/ImmunizedATT.R")

# Data simulation files
source("functions/ClassicDataSim.R") 


### MC XP FUNCTION
MonteCarloSimu <- function(Ry=.8,Rd=.2, n=250,p=300,R=500){
  Results <- matrix(ncol=3, nrow=R)
  AsySD <- matrix(ncol=3, nrow=R)
  t_start <- Sys.time()
  pb <- txtProgressBar(style = 3)
  
  for(r in 1:R){
    ### 1. Generate data
    data <- DataSim(n=n,p=p,Ry=Ry,Rd=Rd,TreatHeter=F)
    X <- data$X
    y <- data$y
    d <- data$d
    
    ### 2. Calibration part
    CAL <- CalibrationLasso(d,X,c=.7,maxIterPen=1e4,PostLasso=T,trace=F)
    
    ### 3. Computes the orthogonality parameter, using method WLS Lasso
    ORT_WLS_L <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="WLSLasso",
                                  c=1.1, nopenset=c(1), RescaleY=F,
                                  maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,PostLasso=F,trace=F)
    
    ### 4. BCH (2014) Estimate
    BCH <- BCHDoubleSelec(y,d,X,cd=.95,cy=1.1,
                          nopenset=c(1),RescaleY=F,
                          maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,trace=F)
    
    
    ### 6. Third step: ATT estimation
    Results[r,] <- c(ImmunizedATT(y,d,X,CAL$betaLasso, Immunity=F)$theta,
                     ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$theta,
                     BCH$theta)
    
    AsySD[r,] <- c(ImmunizedATT(y,d,X,CAL$betaLasso, Immunity=F)$sigma,
                   ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$sigma,
                   BCH$sigma)
    
    setTxtProgressBar(pb, r/R)
    
  }
  
  close(pb)
  print(Sys.time()-t_start)
  return(list(Results=Results, AsySD=AsySD))
}

### RUN SIMULATIONS
i <- 0
z_BEAST <- matrix(,ncol=10,nrow=10)
z_NAIF <- matrix(,ncol=10,nrow=10)
z_BCH <- matrix(,ncol=10,nrow=10)

for(Rd in seq(0,.9,by=.1)){
  i <- i+1
  j <- 0
  print(paste("Rd:",Rd))
  for(Ry in seq(0,.9,by=.1)){
    j <- j+1
    print(paste("Ry:",Ry))
    MCX <- MonteCarloSimu(Ry=Ry,Rd=Rd, n=250,p=300,R=500)
    t <- MCX$Results/MCX$AsySD
    z_NAIF[i,j] <- sum(ifelse(abs(t[,1]) > qnorm(.975),1,0))/500
    z_BEAST[i,j] <- sum(ifelse(abs(t[,2]) > qnorm(.975),1,0))/500
    z_BCH[i,j] <- sum(ifelse(abs(t[,3]) > qnorm(.975),1,0))/500
    print(paste("Naif type I error:",z_NAIF[i,j]))
    print(paste("BEAST type I error:",z_BEAST[i,j]))
    print(paste("BCH type I error:",z_BEAST[i,j]))
  }
}

### Interval coverage plot
x=seq(0,.9,by=.1)

op <- par(bg = "white")
pdf("plots/Confidence_covering_shape.pdf")
persp(x,x,z_BEAST,
      theta=-60,
      phi=10,
      expand=0.5,
      ticktype="detailed",
      col = "steelblue",
      main="BEAST type I error",
      xlim=c(0,1),ylim=c(0,1),zlim=c(0,1),
      xlab="R.2 Selection",
      ylab="R.2 Outcome",
      zlab="Type I error"
)
dev.off()