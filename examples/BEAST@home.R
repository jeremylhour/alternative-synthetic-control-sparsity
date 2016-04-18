### BEAST @ home
### Jeremy L Hour
### 27 mars 2016

### Set working directory
setwd("/Users/jeremylhour/Documents/R/BEAST")

### Load user-defined function
source("functions/CalibrationLasso.R")
source("functions/OrthogonalityReg.R")
source("functions/LogitLasso.R")
source("functions/BCHDoubleSelec.R")
source("functions/VanillaCalibrationLasso.R")
source("functions/VanillaOrthogonalityReg.R")
source("functions/ImmunizedATT.R")

### Load packages
library("AER")
library("lattice")
library("ggplot2")
library("RColorBrewer")

### More Guns, Less Crime example
data("Guns")

head(Guns)

myPalette <- colorRampPalette(rev(brewer.pal(11, "Spectral")), space="Lab")

ggplot(Guns,aes(x = year, y = state, fill = law)) +
          geom_tile() +
          scale_x_discrete(expand = c(0, 0)) +
          scale_y_discrete(expand = c(0, 0)) +
          coord_equal() +
          theme_bw()

fm4 <- lm(log(violent) ~ law + prisoners + density + income +
            population + afam + cauc + male + state + year, data = Guns)
summary(fm4)

### BEAST Model
d <- as.numeric(Guns$law=="yes")
n <- length(d)
X <- Guns[,c("prisoners","density","income","population","afam","cauc","male")]
time.dum <- model.matrix(~ year -1, data=Guns)
time.dum <- time.dum[,-1]
X <- cbind(rep(1,n),X,time.dum)
X <- as.matrix(X)
y <- log(Guns$violent)

CAL <- CalibrationLasso(d,X,c=1.1,maxIterPen=5e1,PostLasso=T,trace=T,maxIter=1e6)
W <- (1-d) * exp(X%*%CAL$betaLasso)/sum(d)
BC <- t(d/sum(d) - W)%*%X

### 3. Computes the orthogonality parameter, using method WLS Lasso
ORT_WLS_L <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="WLSLasso",
                              c=5, nopenset=c(1), RescaleY=F,
                              maxIterPen=1e4,maxIterLasso=1e6,tolLasso=1e-6,PostLasso=F,trace=T)

ORT_WLS_PL <- OrthogonalityReg(y,d,X,CAL$betaPL,method="WLSLasso",
                               c=5, nopenset=c(1), RescaleY=F,
                               maxIterPen=1e4,maxIterLasso=1e6,tolLasso=1e-6,PostLasso=T,trace=T)

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
Results <- c(ImmunizedATT(y,d,X,CAL$betaLasso, Immunity=F)$theta,
                 ImmunizedATT(y,d,X,CAL$betaPL, Immunity=F)$theta,
                 ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$theta,
                 ImmunizedATT(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL, Immunity=T)$theta,
                 ImmunizedATT(y,d,X,LOGIT$betaLasso, Immunity=F)$theta,
                 ImmunizedATT(y,d,X,LOGIT$betaPL, Immunity=F)$theta,
                 BCH$theta,
                 ImmunizedATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso, Immunity=T)$theta,
                 ImmunizedATT(y,d,X,LOGIT$betaPL, FARRELL$muPL, Immunity=T)$theta)

AsySD <- c(ImmunizedATT(y,d,X,CAL$betaLasso, Immunity=F)$sigma,
               ImmunizedATT(y,d,X,CAL$betaPL, Immunity=F)$sigma,
               ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$sigma,
               ImmunizedATT(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL, Immunity=T)$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaLasso, Immunity=F)$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaPL, Immunity=F)$sigma,
               BCH$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso, Immunity=T)$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaPL, FARRELL$muPL, Immunity=T)$sigma)

