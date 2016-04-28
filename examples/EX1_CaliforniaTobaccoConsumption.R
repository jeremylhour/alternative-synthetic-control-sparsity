### Synthetic control example
### Jeremy L Hour and Marianne Blehaut
### 20 janvier 2016
### Last edited: 7 avril 2016

rm(list=ls())

### Set working directory
setwd("R:/Simulations/BEAST")
#setwd("Z:/Simulations/BEAST")

### Load user-defined function
source("functions/CalibrationLasso.R")
source("functions/OrthogonalityReg.R")
source("functions/ImmunizedATT.R")

### Load packages
library("Synth")

####################
####################
####################
### TOBACCO CONS. ##
####################
####################
####################

data <- data.frame(t(read.table("R:/Simulations/MLAB_data.txt")))
#data <- data.frame(t(read.table("Z:/Simulations/MLAB_data.txt")))

Names <- c("State_ID","Income","Retail price", "Young", "BeerCons","Smoking1988", "Smoking1980","Smoking1975",
           mapply(function(x) paste("SmokingCons",x,sep=""),1970:2000))
colnames(data) <- Names
data[,"Treated"] <- as.numeric(data[,"State_ID"]==3) #California is state with ID=3

CaliSmoke <- ts(unlist(data[data[,"Treated"]==1, mapply(function(x) paste("SmokingCons",x,sep=""),1970:2000)]),
                start=c(1970), freq=1)


### Figure 1: California cigarettes consumption
plot(CaliSmoke, col="steelblue", lwd=2,
     xlab="", ylab="Packs per capita",
     ylim=c(0,150),
     main="California cigarette pack consumption")
abline(v=1988, col="firebrick")

### Treatment effect
X <- cbind(rep(1,nrow(data)),
           data[,c("Income","Retail price", "Young", "BeerCons"
                  , "SmokingCons1970"
                  , "SmokingCons1971", "SmokingCons1972", "SmokingCons1973", "SmokingCons1974"
                  , "SmokingCons1975"
                  #, "SmokingCons1976", "SmokingCons1977", "SmokingCons1978", "SmokingCons1979"
                  , "SmokingCons1980"
                  #, "SmokingCons1981", "SmokingCons1982", "SmokingCons1983", "SmokingCons1984"
                  #, "SmokingCons1985", "SmokingCons1986", "SmokingCons1987"
                  , "SmokingCons1988"
                  )])
d <- data[,"Treated"]

### 2. Calibration part

### Setting

d <- as.matrix(d)
X <- as.matrix(X)

n <- nrow(X)
p <- ncol(X)

### First step: Lasso for Calibration
CAL <- CalibrationLasso(d,X,c=0.02,maxIterPen=1e4,trace=T,PostLasso=T)
W <- (1-d)*exp(X%*%CAL$betaLasso)/sum(d)
sum(W)

### Compute ATT
ATT <- vector(length=length(1970:2000))
Immunized <- vector(length=length(1970:2000))
Immunizedsd <- vector(length=length(1970:2000))
NaivePlugIn <- vector(length=length(1970:2000))
Naive <- vector(length=length(1970:2000))

i <- 0
for(t in 1970:2000){
  i=i+1
  varname <- paste("SmokingCons",t,sep="")
  y <- data[,varname]
  
  # Estimation
  ORT <- OrthogonalityReg(y,d,X,beta=CAL$betaLasso, method="WLSLasso",c=.45,
                      maxIterPen=1000,nopenset=1,RescaleY=F,
                      maxIterLasso=10e6,PostLasso=F,trace=F)
  print(paste("year ",t," : "))
  print(ORT$SHat)
  
  #Computing the quantities of interest
  ImmunATT <- ImmunizedATT(y,d,X,CAL$betaLasso,ORT$muLasso, Immunity=T)
  ATT[i] <- ImmunATT$theta
  Immunizedsd[i] <- ImmunATT$sigma
  Immunized[i] <- y[d==1] - ATT[i]
  NaivePlugIn[i] <- y[d==1] - ImmunizedATT(y,d,X,CAL$betaLasso,rep(0,p), Immunity=F)$theta
  Naive[i] <- mean(y[d==0])
}

ATT <- ts(ATT, start=c(1970), freq=1)

ATTdata <- ts(cbind(ATT-1.96*Immunizedsd,ATT,ATT+1.96*Immunizedsd),start=c(1970), freq=1)

### Figure 2: Average Treatment Effect
pdf("plots/Proposition99TreatmentEffect.pdf", width=10, height=6)
plot(ATTdata, plot.type="single",
     col=c("firebrick","firebrick","firebrick"), lwd=c(1,2,1),
     lty=c(6,1,6),xlab="", ylab="Cigarette consumption (Packs per capita)",
     ylim=c(-50,10))
abline(h=0,
       lty=2,col="grey")
lim <- par("usr")
rect(1988, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
legend(1970,-35,
       legend=c("Estimate", ".95 confidence bands"),
       col=c("firebrick","firebrick"), lwd=c(2,1),
       lty=c(1,6))
dev.off()

### Figure 3
plotdata <- ts(cbind(CaliSmoke,Immunized,NaivePlugIn,Naive),start=c(1970), freq=1)

pdf("plots/CaliforniaImmunized.pdf", width=10, height=10)
plot(plotdata, plot.type="single",
     col=c("steelblue","firebrick","forestgreen","darkorchid"), lwd=2,
     lty=c(1,6,6,6),xlab="", ylab="Cigarette consumption (Packs per capita)",
     ylim=c(35,150))
lim <- par("usr")
rect(1988, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
legend(1971,80,
       legend=c("Real California", "Immunized Lasso", "Naive Plug-In Lasso", "Naive Mean"),
       col=c("steelblue","firebrick","forestgreen","darkorchid"), lwd=2,
       lty=c(1,6,6,6))
dev.off()

### 3. Comparison with Synthetic control

# Pre-treatment outcome variables to optimize the fit
Z = NULL
for(t in c(1970:1975, 1980, 1988)){
  varname <- paste("SmokingCons",t,sep="")
  Z <- cbind(Z,data[,varname])
}

# Find the weights...
ADH_SC <- synth(data.prep.obj = NULL,
      X1 = as.matrix(X[d==1,-1]), X0 = t(as.matrix(X[d==0,-1])),
      Z1 = as.matrix(Z[d==1,]), Z0 = t(as.matrix(Z[d==0,])), 
      custom.v = NULL,
      optimxmethod = "Nelder-Mead",
      genoud = FALSE, quadopt = "ipop",
      Margin.ipop = 5e-04,
      Sigf.ipop = 5,
      Bound.ipop = 10,
      verbose = FALSE)

NormalSynth <- synth(data.prep.obj = NULL,
                X1 = as.matrix(X[d==1,-1]), X0 = t(as.matrix(X[d==0,-1])),
                Z1 = as.matrix(X[d==1,-1]), Z0 = t(as.matrix(X[d==0,-1])), 
                custom.v = rep(1,p-1),
                optimxmethod = "Nelder-Mead",
                genoud = FALSE, quadopt = "ipop",
                Margin.ipop = 5e-04,
                Sigf.ipop = 5,
                Bound.ipop = 10,
                verbose = FALSE)

# Construct counterfactual...
ADH <- vector(length=length(1970:2000))
SC <- vector(length=length(1970:2000))
i <- 0
for(t in 1970:2000){
  i=i+1
  varname <- paste("SmokingCons",t,sep="")
  y <- data[,varname]
  ADH[i] <- t(y[d==0])%*%ADH_SC$solution.w
  SC[i] <- t(y[d==0])%*%NormalSynth$solution.w
}

### Figure 4: Synthetic Control plot
plotdata <- ts(cbind(CaliSmoke,Immunized,ADH,SC),start=c(1970), freq=1)

pdf("plots/CaliforniaImmunizedvADH.pdf", width=10, height=10)
plot(plotdata, plot.type="single",
     col=c("steelblue","firebrick","forestgreen","darkgoldenrod1"), lwd=2,
     lty=c(1,6,6,6),xlab="", ylab="Cigarette consumption (Packs per capita)",
     ylim=c(35,150))
lim <- par("usr")
rect(1988, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
legend(1971,80,
       legend=c("Real California", "Immunized, Lasso","Synthetic California (weighted)", "Synthetic California (unweighted)"),
       col=c("steelblue","firebrick","forestgreen","darkgoldenrod1"), lwd=2,
       lty=c(1,6,6,6))
dev.off()


### Balance check
# Checking covariate balancing
# Lasso:
W_Lasso <- (1-d) * exp(X%*%CAL$betaLasso)/sum(d)

# Post-Lasso : 
W_PL <- (1-d) * exp(X%*%CAL$betaPL)/sum(d)

# Synthetic Control
W_SC <- c(ADH_SC$solution.w,0)

BalanceCheck <- data.frame(Method=character(),
                           VarID=double(),
                           mDiff=double(),
                           stringsAsFactors=FALSE)
p <- ncol(X)
VarNames <- c("Constant","Income","Retail price", "Young", "BeerCons", "Smoking1988", "Smoking1980","Smoking1975",
              "SmokingCons1970", "SmokingCons1971", "SmokingCons1972", "SmokingCons1973", "SmokingCons1974")

BalanceCheck[1:p,"Method"] <- rep("Original",p)
BalanceCheck[1:p,"VarID"] <- VarNames 
BalanceCheck[1:p,"mDiff"] <- c(t(d/sum(d) - (1-d)/sum(1-d))%*%X)

BalanceCheck[(p+1):(2*p),"Method"] <- rep("Lasso",p)
BalanceCheck[(p+1):(2*p),"VarID"] <- VarNames 
BalanceCheck[(p+1):(2*p),"mDiff"] <- c(t(d/sum(d) - W_Lasso)%*%X)

BalanceCheck[(2*p+1):(3*p),"Method"] <- rep("PostLasso",p)
BalanceCheck[(2*p+1):(3*p),"VarID"] <- VarNames 
BalanceCheck[(2*p+1):(3*p),"mDiff"] <- c(t(d/sum(d) - W_PL)%*%X)

BalanceCheck[(3*p+1):(4*p),"Method"] <- rep("SyntheticControl",p)
BalanceCheck[(3*p+1):(4*p),"VarID"] <- VarNames 
BalanceCheck[(3*p+1):(4*p),"mDiff"] <- c(t(d/sum(d) - W_SC)%*%X)

BCheck <- subset(BalanceCheck,VarID%in%c("Income","Retail price", "Young", "BeerCons"))
M <- max(abs(BCheck$mDiff))
lb <- 1.1*M
ub <- -1.1*M

pdf("plots/CaliforniaBalanceCheck_X_1980.pdf", width=12, height=12)
ggplot(BCheck, aes(x = mDiff, y = factor(Method, level=c("PostLasso","Lasso","SyntheticControl","Original")))) + 
  geom_point(aes(colour=factor(Method)), size=4) +
  facet_wrap(~VarID, ncol=1) +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none") +
  scale_x_continuous(limits=c(ub,lb),name="Difference in mean, Treated-Control") + 
  ylab("Balancing Method") +
  ggtitle("Balancing check") 
dev.off()

BCheck <- subset(BalanceCheck,VarID%in%c("Smoking1988", "Smoking1980","Smoking1975"))
M <- max(abs(BCheck$mDiff))
lb <- 1.1*M
ub <- -1.1*M

pdf("plots/CaliforniaBalanceCheck_Y_1980.pdf", width=12, height=12)
ggplot(BCheck, aes(x = mDiff, y = factor(Method, level=c("PostLasso","Lasso","SyntheticControl","Original")))) + 
  geom_point(aes(colour=factor(Method)), size=4) +
  facet_wrap(~VarID, ncol=1) +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none") +
  scale_x_continuous(limits=c(ub,lb),name="Difference in mean, Treated-Control") + 
  ylab("Balancing Method") +
  ggtitle("Balancing check") 
dev.off()




#####################
#####################
#####################
### END EXAMPLE #####
#####################
#####################
#####################

### Method 0. Synthetic Control: Original ADH 2010 method
delta <- .0001
v <- rep(1,p) # Weights for importance of each calibration equation
H <- X[d==0,] %*% t(X[d==0,]) + delta * diag(n-1)
f <- c(X[d==1,] %*% t(X[d==0,]))
A <- t(rbind(rep(1,n-1),
           diag(n-1),
           -diag(n-1)))
bc <- c(1,rep(0,n-1),rep(-1,n-1)) 

library("quadprog")
SCsol <- solve.QP(H,f,Amat=A,bvec=bc, meq=1)

SC_CF <- vector(length=length(1970:2000))
i <- 0
for(t in 1970:2000){
  i=i+1
  varname <- paste("SmokingCons",t,sep="")
  y <- data[,varname]
  SC_CF[i] <- t(SCsol$solution)%*%y[d==0]
}

### Figure 4: Synthetic Control plot
plotdata <- ts(cbind(CaliSmoke,SC_CF),start=c(1970), freq=1)

plot(plotdata, plot.type="single",
     col=c("steelblue","darkorchid"), lwd=2,
     lty=c(1,6),xlab="", ylab="Cigarette consumption (Packs per capita)",
     ylim=c(35,150))
lim <- par("usr")
rect(1988, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
legend(1971,80,
       legend=c("Real California", "Ridge-like Regularized Synthetic"),
       col=c("steelblue","darkorchid"), lwd=2,
       lty=c(1,6))



### Method 1. Synthetic Control with Ridge-like regularisation, no constraint
delta <- 1
W_SC <- solve(X[d==0,] %*% t(X[d==0,]) + delta * diag(n-1)) %*% (X[d==0,] %*% as.vector(X[d==1,]))

RidgeCF <- vector(length=length(1970:2000))
i <- 0
for(t in 1970:2000){
  i=i+1
  varname <- paste("SmokingCons",t,sep="")
  y <- data[,varname]
  RidgeCF[i] <- t(W_SC)%*%y[d==0]
}

### Figure 5: Ridge-like regularisation counterfactual
plotdata <- ts(cbind(CaliSmoke,RidgeCF),start=c(1970), freq=1)

plot(plotdata, plot.type="single",
     col=c("steelblue","darkorchid"), lwd=2,
     lty=c(1,6),xlab="", ylab="Cigarette consumption (Packs per capita)",
     ylim=c(35,150))
lim <- par("usr")
rect(1988, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
legend(1971,80,
       legend=c("Real California", "Ridge-like Regularized Synthetic"),
       col=c("steelblue","darkorchid"), lwd=2,
       lty=c(1,6))



### Alternative method 1
### Regression weight counterfactual
### See Abadie, Diamond and Heinmueller (2015)
delta <- 0
W_reg <- as.vector(X[d==1,] %*% solve(t(X[d==0,]) %*% X[d==0,] + delta * diag(p)) %*% t(X[d==0,]))

### Compute Counterfactual
Counterfactual_reg <- vector(length=length(1970:2000))
i <- 0
for(t in 1970:2000){
  i=i+1
  varname <- paste("SmokingCons",t,sep="")
  y <- data[,varname]
  Counterfactual_reg[i] <- W_reg%*%y[d==0]
}

### Figure 6: regression counterfactual
plotdata <- ts(cbind(CaliSmoke,Counterfactual_reg),start=c(1970), freq=1)

pdf("plots/CounterfactualRegression.pdf", width=10, height=10)
plot(plotdata, plot.type="single",
     col=c("steelblue","darkorchid"), lwd=2,
     lty=c(1,6),xlab="", ylab="Cigarette consumption (Packs per capita)",
     ylim=c(35,150))
lim <- par("usr")
rect(1988, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
legend(1971,80,
       legend=c("Real California", "Regression Counterfactual"),
       col=c("steelblue","darkorchid"), lwd=2,
       lty=c(1,6))
dev.off()

### Alternative method 2
### Compute the weight given to each caracteristics by year
### Make the change in the weights not too wigly
Tfin <- length(1970:2000)
bigX <- diag(Tfin) %x% X[d==0,]
bigY <- NULL

for(t in 1970:2000){
  varname <- paste("SmokingCons",t,sep="")
  y <- data[,varname]
  bigY <- c(bigY,y[d==0])
}

### Creation de la matrice de double-difference
ADiff <- matrix(0,nrow=Tfin-2,ncol=Tfin)
for(t in 1:(Tfin-2)){
  ADiff[t,t:(t+2)] <- c(1,-2,1)
}

tildeA <- ADiff %x% t(rep(1,p))

h <- 27
betaSS <- solve(t(bigX)%*%bigX + h*t(tildeA)%*%tildeA) %*% t(bigX) %*% bigY

CF_SS <- t(diag(Tfin) %x% X[d==1,]) %*% betaSS

### Figure 5: Smoothing spline counterfactual
plotdata <- ts(cbind(CaliSmoke,CF_SS),start=c(1970), freq=1)

plot(plotdata, plot.type="single",
     col=c("steelblue","darkorchid"), lwd=2,
     lty=c(1,6),xlab="", ylab="Cigarette consumption (Packs per capita)",
     ylim=c(35,150))
lim <- par("usr")
rect(1988, lim[3], lim[2], lim[4], col = rgb(0.5,0.5,0.5,1/4))
axis(1) ## add axes back
axis(2)
box() 
legend(1971,80,
       legend=c("Real California", "Smoothing spline Counterfactual"),
       col=c("steelblue","darkorchid"), lwd=2,
       lty=c(1,6))


tildeA %*%betaSS