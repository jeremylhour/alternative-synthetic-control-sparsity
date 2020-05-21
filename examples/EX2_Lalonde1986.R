####################
####################
####################
### LALONDE 1986 ###
####################
####################
####################

### Set working directory
setwd("W:/1A_These/A. Research/beast_git/BEAST")
rm(list=ls())

set.seed(9081993)


### 0. Settings
### Load packages
library("ggplot2")
library("lbfgs")
library("causalsens")

### Load user-defined functions
source("functions/CalibrationLasso.R")
source("functions/OrthogonalityReg.R")
source("functions/Compute_ATT.R")
source("functions/LogitLasso.R")
source("functions/BCHDoubleSelec.R")
source("functions/LassoFISTA.R")

### min-max scale
mMscale <- function(X){
  X <- as.matrix(X)
  mins <- apply(X,2,min)
  maxs <- apply(X,2,max)
  return(scale(X, center=mins, scale=maxs-mins))
}

# LOAD AND MANAGE DATA
library("causalsens")
data(lalonde.psid)

d <- lalonde.psid[,"treat"]
y <- lalonde.psid[,"re78"]

### Transformations of the original dataset to have high-dimension data
X.Main <- data.frame("Cons"=1,
                     lalonde.psid[,c("age","education","married","black","hispanic","re74","re75","nodegree")],
                     "NoIncome74"=as.numeric(lalonde.psid[,"re74"]==0),
                     "NoIncome75"=as.numeric(lalonde.psid[,"re75"]==0)
)

X.Main[,c("age","education","re74","re75")] <- mMscale(X.Main[,c("age","education","re74","re75")])

X.ContInterac <- data.frame(
  "AgexMarried"=lalonde.psid[,"age"]*lalonde.psid[,"married"],
  "AgexNodegree"=lalonde.psid[,"age"]*lalonde.psid[,"nodegree"],
  "AgexBlack"=lalonde.psid[,"age"]*lalonde.psid[,"black"],
  "AgexHispanic"=lalonde.psid[,"age"]*lalonde.psid[,"hispanic"],
  "AgexNoIncome74"=lalonde.psid[,"age"]*(lalonde.psid[,"re74"]==0),
  "AgexNoIncome75"=lalonde.psid[,"age"]*(lalonde.psid[,"re75"]==0),
  "EducxMarried"=lalonde.psid[,"education"]*lalonde.psid[,"married"],
  "EducxNodegree"=lalonde.psid[,"education"]*lalonde.psid[,"nodegree"],
  "EducxBlack"=lalonde.psid[,"education"]*lalonde.psid[,"black"],
  "EducxHispanic"=lalonde.psid[,"education"]*lalonde.psid[,"hispanic"],
  "EducxNoIncome74"=lalonde.psid[,"education"]*(lalonde.psid[,"re74"]==0),
  "EducxNoIncome75"=lalonde.psid[,"education"]*(lalonde.psid[,"re75"]==0),
  "Income74xMarried"=lalonde.psid[,"re74"]*lalonde.psid[,"married"],
  "Income74xNoDegree"=lalonde.psid[,"re74"]*lalonde.psid[,"nodegree"],
  "Income74xBlack"=lalonde.psid[,"re74"]*lalonde.psid[,"black"],
  "Income74xHispanic"=lalonde.psid[,"re74"]*lalonde.psid[,"hispanic"],
  "Income74xNoIncome75"=lalonde.psid[,"re74"]*(lalonde.psid[,"re75"]==0),
  "Income75xMarried"=lalonde.psid[,"re75"]*lalonde.psid[,"married"],
  "Income75xNodegree"=lalonde.psid[,"re75"]*lalonde.psid[,"nodegree"],
  "Income75xBlack"=lalonde.psid[,"re75"]*lalonde.psid[,"black"],
  "Income75xHispanic"=lalonde.psid[,"re75"]*lalonde.psid[,"hispanic"],
  "Income75xNoIncome74"=lalonde.psid[,"re75"]*(lalonde.psid[,"re74"]==0)
)

X.ContInterac <- mMscale(X.ContInterac)

X.DumInterac <- data.frame(
  "MarriedxNodegree"=lalonde.psid[,"married"]*lalonde.psid[,"nodegree"],
  "MarriedxBlack"=lalonde.psid[,"married"]*lalonde.psid[,"black"],
  "MarriedxHispanic"=lalonde.psid[,"married"]*lalonde.psid[,"hispanic"],
  "MarriedxNoIncome74"=lalonde.psid[,"married"]*(lalonde.psid[,"re74"]==0),
  "MarriedxNoIncome75"=lalonde.psid[,"married"]*(lalonde.psid[,"re75"]==0),
  "NodegreexBlack"=lalonde.psid[,"nodegree"]*lalonde.psid[,"black"],
  "NodegreexHispanic"=lalonde.psid[,"nodegree"]*lalonde.psid[,"hispanic"],
  "NodegreexNoIncome74"=lalonde.psid[,"nodegree"]*(lalonde.psid[,"re74"]==0),
  "NodegreexNoIncome75"=lalonde.psid[,"nodegree"]*(lalonde.psid[,"re75"]==0),
  "BlackxNoIncome74"=lalonde.psid[,"black"]*(lalonde.psid[,"re74"]==0),
  "BlackxNoIncome75"=lalonde.psid[,"black"]*(lalonde.psid[,"re75"]==0),
  "HispanicxNoIncome74"=lalonde.psid[,"hispanic"]*(lalonde.psid[,"re74"]==0),
  "HispanicxNoIncome75"=lalonde.psid[,"hispanic"]*(lalonde.psid[,"re75"]==0),
  "NoIncome74x75"=(lalonde.psid[,"re74"]==0)*(lalonde.psid[,"re75"]==0)
)

X.Poly <- poly(as.matrix(lalonde.psid[,c("age","education","re74","re75")]),degree=5)
X.Poly <- mMscale(X.Poly)
X.Age <- model.matrix(~as.factor(lalonde.psid[,"age"]) - 1)[,-1]

X <- cbind(X.Main,X.ContInterac,X.DumInterac,X.Poly)
X <- as.matrix(X)


### 2. Balancing
CAL <- CalibrationLasso(d,X,c=1.1,maxIterPen=1e4,PostLasso=T,trace=T)

### 3. Computes the orthogonality parameter
ORT_WLS <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="WLSLasso",
                            c=1.1, nopenset=c(1), RescaleY=F,
                            maxIterPen=1e4,maxIterLasso=1e6,tolLasso=1e-6,PostLasso=F,trace=T)

ORT_WLS_PL <- OrthogonalityReg(y,d,X,CAL$betaPL,method="WLSLasso",
                               c=1.1, nopenset=c(1), RescaleY=F,
                               maxIterPen=1e4,maxIterLasso=1e6,tolLasso=1e-6,PostLasso=T,trace=T)


### 5. Saving results
Results <- data.frame(Estimator=character(),
                      ATT=double(),
                      asymptoticsd=double(),
                      LB95=double(),
                      UB95=double(),
                      PropScore=integer(),
                      Outcome=integer(),
                      stringsAsFactors=FALSE)
### 1. Naive Lasso
Results[1,"Estimator"] <- c("Naive Plug-In -- Lasso")
Results[1,"ATT"] <- Compute_ATT(y,d,X,CAL$betaLasso)$theta
Results[1,"asymptoticsd"] <-  Compute_ATT(y,d,X,CAL$betaLasso)$sigma
Results[1,"PropScore"] <- length(CAL$SHat)
Results[1,"Outcome"] <- 0

### 2. Immunized Lasso
Results[2,"Estimator"] <- c("Immunized -- Lasso")
Results[2,"ATT"] <- Compute_ATT(y,d,X,CAL$betaLasso,ORT_WLS$muLasso)$theta
Results[2,"asymptoticsd"] <- Compute_ATT(y,d,X,CAL$betaLasso,ORT_WLS$muLasso)$sigma
Results[2,"PropScore"] <- length(CAL$SHat)
Results[2,"Outcome"] <- length(ORT_WLS$SHat)

### 3. Naive Post-Lasso
Results[3,"Estimator"] <- c("Naive Plug-In -- Post-Lasso")
Results[3,"ATT"] <- Compute_ATT(y,d,X,CAL$betaPL)$theta
Results[3,"asymptoticsd"] <- Compute_ATT(y,d,X,CAL$betaPL)$sigma
Results[3,"PropScore"] <- length(CAL$SHat)
Results[3,"Outcome"] <- 0

### 4. Immunized Post-Lasso
Results[4,"Estimator"] <- c("Immunized -- Post-Lasso")
Results[4,"ATT"] <- Compute_ATT(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL)$theta
Results[4,"asymptoticsd"] <- Compute_ATT(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL)$sigma
Results[4,"PropScore"] <- length(CAL$SHat)
Results[4,"Outcome"] <- length(ORT_WLS_PL$SHat)


### 6. Other competitors

### Logit Lasso estimate
LOGIT <- LogitLasso(d,X,c=1.1,
                    maxIterPen=1e5,PostLasso=T,trace=T)

### Linear Reg for Farrell (2015)
FARRELL <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="LinearOutcome",
                            c=1.1*sd(y), nopenset=c(1), RescaleY=T,
                            maxIterPen=1e4,maxIterLasso=1e6,tolLasso=1e-6,PostLasso=T,trace=T)

### Pre-selecting the same covariates
### Not faire to compare with our procedure
FARRELL_original <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="LinearOutcome",
                            c=1.1*sd(y), nopenset=c(1,3,7,8), RescaleY=T,
                            maxIterPen=1e4,maxIterLasso=1e6,tolLasso=1e-6,PostLasso=T,trace=T)


### Save Farrell (2015)
Results[5,"Estimator"] <- c("Farrell (2015) -- Lasso")
Results[5,"ATT"] <- Compute_ATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso)$theta
Results[5,"asymptoticsd"] <-  Compute_ATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso)$sigma
Results[5,"PropScore"] <- length(LOGIT$SHat)
Results[5,"Outcome"] <- length(FARRELL$SHat)

Results[6,"Estimator"] <- c("Farrell (2015) -- Post-Lasso")
Results[6,"ATT"] <- Compute_ATT(y,d,X,LOGIT$betaPL, FARRELL$muPL)$theta
Results[6,"asymptoticsd"] <- Compute_ATT(y,d,X,LOGIT$betaPL, FARRELL$muPL)$sigma
Results[6,"PropScore"] <- length(LOGIT$SHat)
Results[6,"Outcome"] <- length(FARRELL$SHat)

### 8. Inverse propensity weighting
Results[9,"Estimator"] <- c("IPW -- Lasso")
Results[9,"ATT"] <- Compute_ATT(y,d,X,LOGIT$betaLasso)$theta
Results[9,"asymptoticsd"] <- Compute_ATT(y,d,X,LOGIT$betaLasso)$sigma
Results[9,"PropScore"] <- length(LOGIT$SHat)
Results[9,"Outcome"] <- 0

Results[10,"Estimator"] <- c("IPW -- Post-Lasso")
Results[10,"ATT"] <- Compute_ATT(y,d,X,LOGIT$betaPL)$theta
Results[10,"asymptoticsd"] <- Compute_ATT(y,d,X,LOGIT$betaPL)$sigma
Results[10,"PropScore"] <- length(LOGIT$SHat)
Results[10,"Outcome"] <- 0

### BCH 2014
BCH <- BCHDoubleSelec(y,d,X,cd=1.1,cy=1.1*sd(y),
                      nopenset=c(1),RescaleY=T,
                      maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,trace=T)

### Save BCH (2014)
Results[7,"Estimator"] <- c("BCH (2014)")
Results[7,"ATT"] <- BCH$theta
Results[7,"asymptoticsd"] <- BCH$sigma
Results[7,"PropScore"] <- length(BCH$SHatd)
Results[7,"Outcome"] <- length(BCH$SHaty)


### 5. Simple linear regression
LinReg <- lm(y ~ d + X)
summary(LinReg)

Results[8,"Estimator"] <- c("OLS, full model")
Results[8,"ATT"] <- coef(LinReg)[2]
Results[8,"asymptoticsd"] <- coef(summary(LinReg))["d","Std. Error"]
Results[8,"PropScore"] <- NA
Results[8,"Outcome"] <- ncol(X)


### Compute confidence interval
Results["LB95"] <- Results[,"ATT"] - qnorm(.975)*Results[,"asymptoticsd"]
Results["UB95"] <- Results[,"ATT"] + qnorm(.975)*Results[,"asymptoticsd"]



### Test of the hdm package
library("hdm")

rY = rlasso(y ~ X)$res
rD = rlasso(d ~ X)$res

print(rD,all=F)

partial.fit.postlasso= lm(rY~rD)
summary(partial.fit.postlasso)$coef["rD",1:2]



################################
################################
################################
# SENSITIVITY TO PENALTY LEVEL #
################################
################################
################################

Yset <- seq(.5,3,by=.2)
Dset <- seq(.3,2,by=.2)

### 5. Saving results
Results <- NULL



for(i in 1:length(Dset)){
  print(i)
	CAL <- CalibrationLasso(d,X,c=Dset[i],maxIterPen=1e4,PostLasso=T,trace=F)

    for(j in 1:length(Yset)){
    print(j)
    ORT_WLS <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="WLSLasso",
                            c=Yset[j], nopenset=c(1), RescaleY=F,
                            maxIterPen=1e4,maxIterLasso=1e6,tolLasso=1e-6,PostLasso=F,trace=F)

### Saving results
### 1. Naive Lasso
NaiveLasso <- c(1,
                ImmunizedATT(y,d,X,CAL$betaLasso,rep(0,p), Immunity=F)$theta,
                ImmunizedATT(y,d,X,CAL$betaLasso,rep(0,p), Immunity=F)$sigma,
                Dset[i],
                Yset[j],
                length(CAL$SHat),
                0)

Results <- rbind(Results, NaiveLasso)

ImmunizedLasso <- c(2,
                    ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS$muLasso, Immunity=T)$theta,
                    ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS$muLasso, Immunity=T)$sigma,
                    Dset[i],
                    Yset[j],
                    length(CAL$SHat),
                    length(ORT_WLS$SHat))

Results <- rbind(Results, ImmunizedLasso)
print(ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS$muLasso, Immunity=T)$theta)

  }
}


res <- data.frame(Results)
names(res) <- c("BEAST","ATT","AsySD","cd","cy","Sd","Sy")
res[,"tstat"] <- res[,"ATT"]/res[,"AsySD"]


library("reshape")
ATT.m <- melt(ATT)
ATT.m[,"Delta.Y"] <- Yset[ATT.m[,1]]
ATT.m[,"Delta.D"] <- Dset[ATT.m[,2]]

pdf("plots/PenaltyHeatMap_BEAST.pdf", height=10, width=10)
p <- ggplot(subset(res, BEAST==2), aes(x=cy,y=cd,fill=ATT)) + 
              geom_tile(aes(fill=ATT), colour="white") +
              scale_fill_gradient(low = "white", high="steelblue",
                                  name="ATT estimate") +
              scale_colour_grey() + theme_bw()

base_size <- 10
p + theme_bw(base_size = base_size) +
  labs(x = "Penalty for Outcome regression",y = "Penalty for Propensity Score", size=2) + 
  theme(legend.position = "right",axis.ticks = element_blank(),
        axis.text.x = element_text(size = base_size * 1.1, hjust = 0, colour = "grey50"),
        axis.text.y = element_text(size = base_size * 1.1, hjust = 0, colour = "grey50"),
        panel.grid.major = element_blank(),panel.grid.minor = element_blank(),
        panel.border = element_blank())
dev.off()



plot(res[res$BEAST==1,"ATT"],res[res$BEAST==2,"ATT"])



### 6. CBPS
library("CBPS")
CBPSfit <- CBPS(treat ~ X,
     data = lalonde.psid, ATT = 1, iterations = 1000,
     standardize = TRUE, method = "over", twostep = TRUE)
summary(CBPSfit)

## Horwitz-Thompson estimate
mean(d*y/CBPSfit$fitted.values)
## Inverse propensity score weighting
sum(d*y/CBPSfit$fitted.values)/sum(d/CBPSfit$fitted.values)

library(MatchIt)
m.out <- matchit(treat ~ fitted(CBPSfit), method = "nearest", data = lalonde.psid,
                 replace = TRUE)
summary(m.out)
m.data <- match.data(m.out)
mean(lalonde.psid[lalonde.psid$treat==1,"re78"])-mean(m.data[,"re78"])