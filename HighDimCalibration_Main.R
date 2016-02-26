### High-dimension Calibration for Treatment Effect
### 28 decembre 2015
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
source("functions/CalibrationLasso.R")
source("functions/OrthogonalityReg.R")
source("functions/LogitLasso.R")
source("functions/BCHDoubleSelec.R")
source("functions/VanillaCalibrationLasso.R")
source("functions/VanillaOrthogonalityReg.R")
source("functions/ImmunizedATT.R")


### MC XP
R <- 1000
Results <- matrix(ncol=9, nrow=R)
AsySD <- matrix(ncol=9, nrow=R)
Convergence <- matrix(ncol=3, nrow=R)
Weights <- Convergence
t_start <- Sys.time()

for(r in 1:R){
### 1. Generate data
print(paste("MC Iteration nb.",r))
  
data <- DataSim(n=200,p=200,Ry=.8,Rd=.2)
X <- data$X
y <- data$y
d <- data$d

### 2. Calibration part
CAL <- CalibrationLasso(d,X,c=.7,maxIterPen=5e1,PostLasso=T,trace=F,maxIter=1e6)
W <- (1-d) * exp(X%*%CAL$betaPL)/sum(d)
BC <- t(d/sum(d) - W)%*%X

### 3. Computes the orthogonality parameter, using method WLS Lasso
ORT_WLS_L <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="WLSLasso",
                              c=2, nopenset=c(1), RescaleY=F,
                              maxIterPen=1e4,maxIterLasso=1e3,tolLasso=1e-6,PostLasso=F,trace=F)

ORT_WLS_PL <- OrthogonalityReg(y,d,X,CAL$betaPL,method="WLSLasso",
                               c=2, nopenset=c(1), RescaleY=F,
                               maxIterPen=1e4,maxIterLasso=1e3,tolLasso=1e-6,PostLasso=T,trace=F)

### 4. Logit Lasso estimate
LOGIT <- LogitLasso(d,X,c=.6,
                    maxIterPen=5e1,PostLasso=T,trace=F)

### 4 bis. Farrell (2015)
FARRELL <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="LinearOutcome",
                              c=2, nopenset=c(1), RescaleY=F,
                              maxIterPen=1e4,maxIterLasso=1e3,tolLasso=1e-6,PostLasso=T,trace=T)

### 5. BCH (2014) Estimate
BCH <- BCHDoubleSelec(y,d,X,cd=.95,cy=1.1,
                       nopenset=c(1),RescaleY=F,
                      maxIterPen=1e4,maxIterLasso=1e3,tolLasso=1e-6,trace=F)

### 6. Third step: ATT estimation
Results[r,] <- c(ImmunizedATT(y,d,X,CAL$betaLasso, Immunity=F)$theta,
                 ImmunizedATT(y,d,X,CAL$betaPL, Immunity=F)$theta,
                 ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$theta,
                 ImmunizedATT(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL, Immunity=T)$theta,
                 ImmunizedATT(y,d,X,LOGIT$betaLasso, Immunity=F)$theta,
                 ImmunizedATT(y,d,X,LOGIT$betaPL, Immunity=F)$theta,
                 BCH$theta,
                 ImmunizedATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso, Immunity=T)$theta,
                 ImmunizedATT(y,d,X,LOGIT$betaPL, FARRELL$muPL, Immunity=T)$theta)

AsySD[r,] <- c(ImmunizedATT(y,d,X,CAL$betaLasso, Immunity=F)$sigma,
               ImmunizedATT(y,d,X,CAL$betaPL, Immunity=F)$sigma,
               ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$sigma,
               ImmunizedATT(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL, Immunity=T)$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaLasso, Immunity=F)$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaPL, Immunity=F)$sigma,
               BCH$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso, Immunity=T)$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaPL, FARRELL$muPL, Immunity=T)$sigma)

Convergence[r,] <- c(CAL$convergence,
                     ORT_WLS_L$convergence,
                     ORT_WLS_PL$convergence)
}

print(Sys.time()-t_start)

# Post-simulation treatment
# Discard all draws which did not converge
valid <- Convergence[,1] == 0 & Convergence[,2]==0  & Convergence[,3] == 0
Results <- Results[valid,]
Results <- Results[!is.na(Results[,1]),]
R <- nrow(Results)

# Draw the charts
id <- c(mapply(function(x) rep(x,R),1:9))
val <- c(Results)
data_res <- data.frame(val = val, model = id)

M <- max(abs(quantile(Results,.01)),abs(quantile(Results,.99)))
lb <- -1.1*M
ub <- 1.1*M


sdBCH <- mean(AsySD[,3])

### Function for plot
get.plot <- function(data,modelS,title="A Title",sdBCH){
  plot_res <- ggplot(subset(data, (model==modelS)), aes(x=val)) + 
    geom_histogram(binwidth = .02, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
    scale_x_continuous(limits=c(lb,ub), name="Treatment effect") +
    ggtitle(title) + 
    stat_function(fun = dnorm, args=list(mean=0, sd=sdBCH), colour="darkorchid3", size=1) +
    theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none")
  
  return(plot_res)
}


pdf("plots/PlugInvsImmunizedLasso.pdf", width=12, height=4)
grid.arrange(get.plot(data_res,1,"Naive Plug-In, Lasso", sdBCH), get.plot(data_res,3,"Immunized, Lasso", sdBCH), ncol=2)
dev.off()

pdf("plots/PlugInvsImmunizedPostLasso.pdf", width=12, height=4)
grid.arrange(get.plot(data_res,2,"Naive Plug-In, Post-Lasso", sdBCH), get.plot(data_res,4,"Immunized, Post-Lasso", sdBCH), ncol=2)
dev.off()

pdf("plots/LogitLasso.pdf", width=12, height=4)
grid.arrange(get.plot(data_res,5,"IPW Logit, Lasso", sdBCH), get.plot(data_res,6,"IPW Logit, Post-Lasso", sdBCH), ncol=2)
dev.off()

pdf("plots/Farrell.pdf", width=12, height=4)
grid.arrange(get.plot(data_res,8,"Farrell, Lasso", sdBCH), get.plot(data_res,9,"Farrell, Post-Lasso", sdBCH), ncol=2)
dev.off()

### Compute bias and RMSE
StatDisplay <- data.frame()
StatDisplay[1:9,"bias"] <- 100*apply(Results,2,mean)
StatDisplay[1:9,"RMSE"]  <- 100*sqrt(apply(Results^2,2,mean))
StatDisplay[1:9,"AsySD"]  <- apply(AsySD,2,mean)
StatDisplay[1:9,"ShapiroTest"]  <- apply(Results,2, function(x) shapiro.test(x)$p.value)
row.names(StatDisplay) <- c("NPILasso","NPIPL","ImmunizedLasso","ImmunizedPL","LogitLasso", "LogitPL", "BCH 2014", "Farrell, Lasso", "Farrell, Post-Lasso")
print(StatDisplay)