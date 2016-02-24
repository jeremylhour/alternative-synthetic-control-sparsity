### High-dimension Calibration for Treatment Effect
### 2 fevrier 2015
### M Blehaut

### Set working directory
setwd("Z:/Simulations/BEAST")
# R pour Jeremy, Z pour Marianne


rm(list=ls())
set.seed(30031987)



### 0. Settings

### Load packages
library("ggplot2")
library("gridExtra")
library("MASS")

### Load user-defined functions
source("functions/DataSim.R") 
source("functions/CalibrationLasso.R")
source("functions/OrthogonalityReg.R")
source("functions/ImmunizedATT.R")


### MC XP
Simu <- function(Rep=1000,N,P,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005){
R <- Rep
Results <- matrix(ncol=3, nrow=R)
Convergence <- Results
Weights <- NULL
B <- matrix(ncol=P+1, nrow=R)
ML <- matrix(ncol=P+1, nrow=R)
MPL <- matrix(ncol=P+1, nrow=R)
t_start <- Sys.time()

for(r in 1:R){
	### 1. Generate data
	print(paste("MC Iteration nb.",r))
  
	data <- DataSim(n=N,p=P,Ry=R2y,Rd=R2d)
	X <- data$X
	y <- data$y
	d <- data$d

	### 2. Calibration part
	CAL <- CalibrationLasso(d,X,c=c1,maxIterPen=5e1,PostLasso=T,trace=F,maxIter=1e6)
	W <- (1-d) * exp(X%*%CAL$betaPL)/sum(d)
	BC <- t(d/sum(d) - W)%*%X

	### 3. Computes the orthogonality parameter, using method WLS Lasso
	ORT_WLS_L <- OrthogonalityReg(y,d,X,beta=CAL$betaLasso,method="WLSLasso",c=c2,
                              maxIterPen=1e5,maxIterLasso=1e6,PostLasso=F,trace=F)
	ORT_WLS_PL <- OrthogonalityReg(y,d,X,beta=CAL$betaPL,method="WLSLasso",c=c2,
                               maxIterPen=1e5,maxIterLasso=1e6,PostLasso=T,trace=F)

	### 5. Third step: ATT estimation
	Results[r,] <- c(ImmunizedATT(y,d,X,CAL$betaLasso,rep(0,p), Immunity=F),
                 ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T),
                 ImmunizedATT(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL, Immunity=T))

	Convergence[r,] <- c(CAL$convergence,
                     ORT_WLS_L$convergence,
                     ORT_WLS_PL$convergence)

	Weights <- cbind(Weights,W)
	B[r,] <- CAL$betaLasso
	ML[r,] <- ORT_WLS_L$muLasso
	MPL[r,] <- ORT_WLS_PL$muPL
	}
tps <- Sys.time()-t_start
print(tps)

# Post-simulation treatment
# Discard all those who did not converge
valid <- Convergence[,1] == 0 & Convergence[,2] == 0 & Convergence[,3] == 0
results <- Results[valid,]
results <- results[!is.na(results[,1]),]
R <- nrow(results)

id <- c(rep(1,R),rep(2,R),rep(3,R))
val <- c(results)
data_res <- data.frame(val = val, model = id)

M <- max(abs(quantile(results,.01)),abs(quantile(results,.99)))
lb <- -1.1*M
ub <- 1.1*M

sdBCH <- sd(results[,2])

plot1 <- ggplot(subset(data_res, (model==1)), aes(x=val)) + 
  geom_histogram(binwidth = .02, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
  scale_x_continuous(limits=c(lb,ub), name="Treatment effect") +
  ggtitle("Naive Plug-In") + 
  stat_function(fun = dnorm, args=list(mean=0, sd=sdBCH), colour="darkorchid3", size=1) +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none") 

plot2 <- ggplot(subset(data_res, (model==2)), aes(x=val)) + 
  geom_histogram(binwidth = .02, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
  scale_x_continuous(limits=c(lb,ub), name="Treatment effect") +
  ggtitle("Immunized, Lasso") + 
  stat_function(fun = dnorm, args=list(mean=0, sd=sdBCH), colour="darkorchid3", size=1) +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none") 

plot3 <- ggplot(subset(data_res, (model==3)), aes(x=val)) + 
  geom_histogram(binwidth = .02, alpha=.5, position='identity',fill="steelblue", aes(y = ..density..)) +
  scale_x_continuous(limits=c(lb,ub), name="Treatment effect") +
  ggtitle("Immunized, Post-Lasso") + 
  stat_function(fun = dnorm, args=list(mean=0, sd=sdBCH), colour="darkorchid3", size=1) +
  theme(plot.title = element_text(lineheight=.8, face="bold"),legend.position="none") 


pdf(paste("plots/Plots_Simus/PlugInvsImmunized_",N,"_",P,".pdf",sep=""), 
width=12, height=4)
grid.arrange(plot1, plot2, plot3, ncol=3)
dev.off()
pdf(paste("plots/Plots_Simus/PlugInvsImmunized_",N,"_",P,"_sple.pdf",sep=""), 
width=8, height=4)
grid.arrange(plot1, plot2, ncol=2)
dev.off()

return(list(Results, Convergence, B, ML, MPL))
}

N50P50 <- Simu(Rep=1000,N=50,P=50,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005)
N50P100 <- Simu(Rep=1000,N=50,P=100,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005)

N100P50 <- Simu(Rep=1000,N=100,P=50,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005)
N100P100 <- Simu(Rep=1000,N=100,P=100,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005)
N100P200 <- Simu(Rep=1000,N=100,P=200,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005)

N200P50 <- Simu(Rep=1000,N=200,P=50,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005)
N200P100 <- Simu(Rep=1000,N=200,P=100,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005)
N200P200 <- Simu(Rep=1000,N=200,P=200,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005)
N200P500 <- Simu(Rep=1000,N=200,P=500,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005)

N500P50 <- Simu(Rep=1000,N=500,P=50,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005)
N500P100 <- Simu(Rep=1000,N=500,P=100,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005)
N500P200 <- Simu(Rep=1000,N=500,P=200,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005)
N500P500 <- Simu(Rep=1000,N=500,P=500,R2y=0.8,R2d=0.2,c1=0.7,c2=0.005)

#save.image("Z:\\Simulations\\R_Code\\simus")

load("Z:\\Simulations\\R_Code\\simus")

res <- matrix(ncol=8, nrow=8)

valid <- N50P50[[2]][,1] == 0 & N50P50[[2]][,2] == 0 & N50P50[[2]][,3] == 0
results <- N50P50[[1]][valid,]
results <- results[!is.na(results[,1]),]
res_pi_50 <- c(mean(results[,1]), sqrt(mean(results[,1]^2)))
res_im_50 <- c(mean(results[,2]), sqrt(mean(results[,2]^2)))

valid <- N50P100[[2]][,1] == 0 & N50P100[[2]][,2] == 0 & N50P100[[2]][,3] == 0
results <- N50P100[[1]][valid,]
results <- results[!is.na(results[,1]),]
res_pi_100 <- c(mean(results[,1]), sqrt(mean(results[,1]^2)))
res_im_100 <- c(mean(results[,2]), sqrt(mean(results[,2]^2)))

res[1,] <- c(res_pi_50, res_pi_100, NA, NA, NA, NA)
res[2,] <- c(res_im_50, res_im_100, NA, NA, NA, NA)

valid <- N100P50[[2]][,1] == 0 & N100P50[[2]][,2] == 0 & N100P50[[2]][,3] == 0
results <- N100P50[[1]][valid,]
results <- results[!is.na(results[,1]),]
res_pi_50 <- c(mean(results[,1]), sqrt(mean(results[,1]^2)))
res_im_50 <- c(mean(results[,2]), sqrt(mean(results[,2]^2)))

valid <- N100P100[[2]][,1] == 0 & N100P100[[2]][,2] == 0 & N100P100[[2]][,3] == 0
results <- N100P100[[1]][valid,]
results <- results[!is.na(results[,1]),]
res_pi_100 <- c(mean(results[,1]), sqrt(mean(results[,1]^2)))
res_im_100 <- c(mean(results[,2]), sqrt(mean(results[,2]^2)))

valid <- N100P200[[2]][,1] == 0 & N100P200[[2]][,2] == 0 & N100P200[[2]][,3] == 0
results <- N100P200[[1]][valid,]
results <- results[!is.na(results[,1]),]
res_pi_200 <- c(mean(results[,1]), sqrt(mean(results[,1]^2)))
res_im_200 <- c(mean(results[,2]), sqrt(mean(results[,2]^2)))

res[3,] <- c(res_pi_50, res_pi_100, res_pi_200, NA, NA)
res[4,] <- c(res_im_50, res_im_100, res_im_200, NA, NA)

valid <- N200P50[[2]][,1] == 0 & N200P50[[2]][,2] == 0 & N200P50[[2]][,3] == 0
results <- N200P50[[1]][valid,]
results <- results[!is.na(results[,1]),]
res_pi_50 <- c(mean(results[,1]), sqrt(mean(results[,1]^2)))
res_im_50 <- c(mean(results[,2]), sqrt(mean(results[,2]^2)))

valid <- N200P100[[2]][,1] == 0 & N200P100[[2]][,2] == 0 & N200P100[[2]][,3] == 0
results <- N200P100[[1]][valid,]
results <- results[!is.na(results[,1]),]
res_pi_100 <- c(mean(results[,1]), sqrt(mean(results[,1]^2)))
res_im_100 <- c(mean(results[,2]), sqrt(mean(results[,2]^2)))

valid <- N200P200[[2]][,1] == 0 & N200P200[[2]][,2] == 0 & N200P200[[2]][,3] == 0
results <- N200P200[[1]][valid,]
results <- results[!is.na(results[,1]),]
res_pi_200 <- c(mean(results[,1]), sqrt(mean(results[,1]^2)))
res_im_200 <- c(mean(results[,2]), sqrt(mean(results[,2]^2)))

valid <- N200P500[[2]][,1] == 0 & N200P500[[2]][,2] == 0 & N200P500[[2]][,3] == 0
results <- N200P500[[1]][valid,]
results <- results[!is.na(results[,1]),]
res_pi_500 <- c(mean(results[,1]), sqrt(mean(results[,1]^2)))
res_im_500 <- c(mean(results[,2]), sqrt(mean(results[,2]^2)))

res[5,] <- c(res_pi_50, res_pi_100, res_pi_200, res_pi_500)
res[6,] <- c(res_im_50, res_im_100, res_im_200, res_pi_500)

valid <- N500P50[[2]][,1] == 0 & N500P50[[2]][,2] == 0 & N500P50[[2]][,3] == 0
results <- N500P50[[1]][valid,]
results <- results[!is.na(results[,1]),]
res_pi_50 <- c(mean(results[,1]), sqrt(mean(results[,1]^2)))
res_im_50 <- c(mean(results[,2]), sqrt(mean(results[,2]^2)))

valid <- N500P100[[2]][,1] == 0 & N500P100[[2]][,2] == 0 & N500P100[[2]][,3] == 0
results <- N500P100[[1]][valid,]
results <- results[!is.na(results[,1]),]
res_pi_100 <- c(mean(results[,1]), sqrt(mean(results[,1]^2)))
res_im_100 <- c(mean(results[,2]), sqrt(mean(results[,2]^2)))

valid <- N500P200[[2]][,1] == 0 & N500P200[[2]][,2] == 0 & N500P200[[2]][,3] == 0
results <- N500P200[[1]][valid,]
results <- results[!is.na(results[,1]),]
res_pi_200 <- c(mean(results[,1]), sqrt(mean(results[,1]^2)))
res_im_200 <- c(mean(results[,2]), sqrt(mean(results[,2]^2)))

valid <- N500P500[[2]][,1] == 0 & N500P500[[2]][,2] == 0 & N500P500[[2]][,3] == 0
results <- N500P500[[1]][valid,]
results <- results[!is.na(results[,1]),]
res_pi_500 <- c(mean(results[,1]), sqrt(mean(results[,1]^2)))
res_im_500 <- c(mean(results[,2]), sqrt(mean(results[,2]^2)))

res[7,] <- c(res_pi_50, res_pi_100, res_pi_200, res_pi_500)
res[8,] <- c(res_im_50, res_im_100, res_im_200, res_pi_500)

res <- format(round(res*100, digits=3), digits=3, trim=T)

res[is.na(as.numeric(res))] <- "-999"
res[as.numeric(res[,1])>0,1] <- paste("~",res[as.numeric(res[,1])>0,1], sep="")
res[as.numeric(res[,3])>0,3] <- paste("~",res[as.numeric(res[,3])>0,3], sep="")
res[as.numeric(res[,5])>0,5] <- paste("~",res[as.numeric(res[,5])>0,5], sep="")
res[as.numeric(res[,7])>0,7] <- paste("~",res[as.numeric(res[,7])>0,7], sep="")
res[res=="-999"] <- NA

write.table(res, file = "table_simu.txt", append = FALSE, quote = FALSE, sep = " & ",
            eol = paste(" \\\\ \n"), na = "--", dec = ".", row.names = F,
            col.names = F)








