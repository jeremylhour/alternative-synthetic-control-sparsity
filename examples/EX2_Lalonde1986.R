####################
####################
####################
### LALONDE 1986 ###
####################
####################
####################

### Set working directory
setwd("R:/Simulations/BEAST") # R pour Jeremy, Z pour Marianne

rm(list=ls())
set.seed(9081993)


### 0. Settings
### Load packages
library("plyr")
library("ggplot2")

### Load user-defined functions
source("functions/CalibrationLasso.R")
source("functions/OrthogonalityReg.R")
source("functions/ImmunizedATT.R")
source("functions/LogitLasso.R")
source("functions/BCHDoubleSelec.R")

# Load data
{ 
  library("causalsens")
data(lalonde.psid)

d <- lalonde.psid[,"treat"]
y <- lalonde.psid[,"re78"]

### Transformations of the original dataset to have high-dimension data
### Rescale continuous variables
X.Main <- data.frame("Cons"=1,
                     lalonde.psid[,c("age","education","married","black","hispanic","re74","re75","nodegree")],
                     "NoIncome74"=as.numeric(lalonde.psid[,"re74"]==0),
                     "NoIncome75"=as.numeric(lalonde.psid[,"re75"]==0)
)

# X.Main[,c("age","education","re74","re75")] <- scale(X.Main[,c("age","education","re74","re75")])

X.ContInterac <- data.frame(
  "AgexMarried"=lalonde.psid[,"age"]*lalonde.psid[,"married"],
  "AgexNodegree"=lalonde.psid[,"age"]*lalonde.psid[,"nodegree"],
  "AgexBlack"=lalonde.psid[,"age"]*lalonde.psid[,"black"],
  "AgexHispanic"=lalonde.psid[,"age"]*lalonde.psid[,"hispanic"],
  "AgexNoIncome74"=lalonde.psid[,"age"]*(lalonde.psid[,"re74"]==0),
  "AgexNoIncome75"=lalonde.psid[,"age"]*(lalonde.psid[,"re75"]==0),
  "EducxMarried"=lalonde.psid[,"education"]*lalonde.psid[,"married"],
  "EdcuxNodegree"=lalonde.psid[,"education"]*lalonde.psid[,"nodegree"],
  "EdcuxBlack"=lalonde.psid[,"education"]*lalonde.psid[,"black"],
  "EdcuxHispanic"=lalonde.psid[,"education"]*lalonde.psid[,"hispanic"],
  "EdcuxNoIncome74"=lalonde.psid[,"education"]*(lalonde.psid[,"re74"]==0),
  "EdcuxNoIncome75"=lalonde.psid[,"education"]*(lalonde.psid[,"re75"]==0),
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

# X.ContInterac <- scale(X.ContInterac)

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
X.Poly <- scale(X.Poly)
# X.Educ <- model.matrix(~as.factor(lalonde.psid[,"education"] - 1))[,-1]
X.Age <- model.matrix(~as.factor(lalonde.psid[,"age"]) - 1)[,-1]

X <- cbind(X.Main,X.ContInterac,X.DumInterac,X.Poly)
X <- as.matrix(X)
}


### 2. Calibration part
CAL <- CalibrationLasso(d,X,c=1.1,maxIterPen=1e4,PostLasso=T,trace=T,maxIter=1e6)

# Checking covariate balancing
# Lasso:
W_Lasso <- (1-d) * exp(X%*%CAL$betaLasso)/sum(d)

# Post-Lasso : 
W_PL <- (1-d) * exp(X%*%CAL$betaPL)/sum(d)

DVar <- 100*(W_PL - W_Lasso)/W_Lasso

plot(W_Lasso,DVar)
abline(h=0,lty=6)

# data <- list(X=X,y=y,d=d,W=W_Lasso)
# save(data,file="Lalondedata_Unscaled.R")

BalanceCheck <- data.frame(Method=character(),
                           VarID=double(),
                           mDiff=double(),
                           stringsAsFactors=FALSE)
p <- ncol(X)

BalanceCheck[1:p,"Method"] <- rep("Original",p)
BalanceCheck[1:p,"VarID"] <- names(cbind(X.Main,X.ContInterac,X.DumInterac,X.Poly))
BalanceCheck[1:p,"mDiff"] <- c(t(d/sum(d) - (1-d)/sum(1-d))%*%X)

BalanceCheck[(p+1):(2*p),"Method"] <- rep("Lasso",p)
BalanceCheck[(p+1):(2*p),"VarID"] <- names(cbind(X.Main,X.ContInterac,X.DumInterac,X.Poly))
BalanceCheck[(p+1):(2*p),"mDiff"] <- c(t(d/sum(d) - W_Lasso)%*%X)

BalanceCheck[(2*p+1):(3*p),"Method"] <- rep("PostLasso",p)
BalanceCheck[(2*p+1):(3*p),"VarID"] <- names(cbind(X.Main,X.ContInterac,X.DumInterac,X.Poly))
BalanceCheck[(2*p+1):(3*p),"mDiff"] <- c(t(d/sum(d) - W_PL)%*%X)

BCheck <- subset(BalanceCheck,VarID%in%c("age","education","married","black","hispanic"))
M <- max(abs(BCheck$mDiff))
lb <- 1.1*M
ub <- -1.1*M

pdf("plots/LalondeBalanceCheck.pdf", width=12, height=12)
ggplot(BCheck, aes(x = mDiff, y = factor(Method, level=c("PostLasso","Lasso","Original")))) + 
  geom_point(aes(colour=factor(Method)), size=4) +
  facet_wrap(~VarID, ncol=1) +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none") +
  scale_x_continuous(limits=c(ub,lb),name="Difference in mean, Treated-Control") + 
  ylab("Balancing Method") +
  ggtitle("Balancing check") 
dev.off()


### 3. Computes the orthogonality parameter
ORT_WLS <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="WLSLasso",
                            c=3*sd(y), nopenset=c(1), RescaleY=T,
                            maxIterPen=1e4,maxIterLasso=1e6,tolLasso=1e-6,PostLasso=F,trace=T)

ORT_WLS_PL <- OrthogonalityReg(y,d,X,CAL$betaPL,method="WLSLasso",
                               c=3*sd(y), nopenset=c(1), RescaleY=T,
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
Results[1,"Estimator"] <- c("Naive Plug-In, Lasso")
Results[1,"ATT"] <- ImmunizedATT(y,d,X,CAL$betaLasso,rep(0,p), Immunity=F)$theta
Results[1,"asymptoticsd"] <-  ImmunizedATT(y,d,X,CAL$betaLasso,rep(0,p), Immunity=F)$sigma
Results[1,"PropScore"] <- length(CAL$SHat)
Results[1,"Outcome"] <- 0

### 2. Immunized Lasso
Results[2,"Estimator"] <- c("Immunized, Lasso")
Results[2,"ATT"] <- ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS$muLasso, Immunity=T)$theta
Results[2,"asymptoticsd"] <- ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS$muLasso, Immunity=T)$sigma
Results[2,"PropScore"] <- length(CAL$SHat)
Results[2,"Outcome"] <- length(ORT_WLS$SHat)

### 3. Naive Post-Lasso
Results[3,"Estimator"] <- c("Naive Plug-In, Post-Lasso")
Results[3,"ATT"] <- ImmunizedATT(y,d,X,CAL$betaPL,rep(0,p), Immunity=F)$theta
Results[3,"asymptoticsd"] <- ImmunizedATT(y,d,X,CAL$betaPL,rep(0,p), Immunity=F)$sigma
Results[3,"PropScore"] <- length(CAL$SHat)
Results[3,"Outcome"] <- 0

### 4. Immunized Post-Lasso
Results[4,"Estimator"] <- c("Immunized, Post-Lasso")
Results[4,"ATT"] <- ImmunizedATT(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL, Immunity=T)$theta
Results[4,"asymptoticsd"] <- ImmunizedATT(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL, Immunity=T)$sigma
Results[4,"PropScore"] <- length(CAL$SHat)
Results[4,"Outcome"] <- length(ORT_WLS_PL$SHat)


### 6. Other competitors

### Logit Lasso estimate
LOGIT <- LogitLasso(d,X,c=1.001,
                    maxIterPen=1e5,PostLasso=T,trace=T,maxIter=1e6)

### Linear Reg for Farrell (2015)
FARRELL <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="LinearOutcome",
                            c=3*sd(y), nopenset=c(1), RescaleY=T,
                            maxIterPen=1e4,maxIterLasso=1e6,tolLasso=1e-6,PostLasso=T,trace=T)

### Save Farrell (2015)
Results[5,"Estimator"] <- c("Farrell (2015), Lasso")
Results[5,"ATT"] <- ImmunizedATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso, Immunity=T)$theta
Results[5,"asymptoticsd"] <- ImmunizedATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso, Immunity=T)$sigma
Results[5,"PropScore"] <- length(LOGIT$SHat)
Results[5,"Outcome"] <- length(FARRELL$SHat)

Results[6,"Estimator"] <- c("Farrell (2015), Post-Lasso")
Results[6,"ATT"] <- ImmunizedATT(y,d,X,LOGIT$betaPL, FARRELL$muPL, Immunity=T)$theta
Results[6,"asymptoticsd"] <- ImmunizedATT(y,d,X,LOGIT$betaPL, FARRELL$muPL, Immunity=T)$sigma
Results[6,"PropScore"] <- length(LOGIT$SHat)
Results[6,"Outcome"] <- length(FARRELL$SHat)

### BCH 2014
BCH <- BCHDoubleSelec(y,d,X,cd=2,cy=1.1,
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


### 8. Inverse propensity weighting
Results[9,"Estimator"] <- c("IPW, Lasso")
Results[9,"ATT"] <- ImmunizedATT(y,d,X,LOGIT$betaLasso, Immunity=F)$theta
Results[9,"asymptoticsd"] <- ImmunizedATT(y,d,X,LOGIT$betaLasso, Immunity=F)$sigma
Results[9,"PropScore"] <- length(LOGIT$SHat)
Results[9,"Outcome"] <- 0

Results[10,"Estimator"] <- c("IPW, Post-Lasso")
Results[10,"ATT"] <- ImmunizedATT(y,d,X,LOGIT$betaPL, Immunity=F)$theta
Results[10,"asymptoticsd"] <- ImmunizedATT(y,d,X,LOGIT$betaPL, Immunity=F)$sigma
Results[10,"PropScore"] <- length(LOGIT$SHat)
Results[10,"Outcome"] <- 0




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