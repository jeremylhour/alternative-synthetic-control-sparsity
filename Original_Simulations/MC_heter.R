### High-dimension Calibration for Treatment Effect
### 2 fevrier 2015
### M Blehaut

### Set working directory
setwd("W:/Bureau/calage_temp/BEAST")


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


### MC XP

Simu_heter <- function(Rep=1000,N,P,R2y=0.8,R2d=0.2,c1=0.7,c2=2){
  
  ## STEP A. SIMULATIONS
  
  R <- Rep
  Results <- matrix(ncol=5, nrow=R)
  AsySD <- matrix(ncol=5, nrow=R)
  Convergence <- matrix(ncol=2, nrow=R)
  t_start <- Sys.time()
  pb <- txtProgressBar(style = 3)
  
  t_start <- Sys.time()
  
  for(r in 1:R){
    print(paste("MC Iteration nb.",r))
    
    ### 1. Generate data
    data <- DataSim(n=N,p=P,Ry=R2y,Rd=R2d,TreatHeter=T)
    X <- data$X
    y <- data$y
    d <- data$d
    
    ### 2. Calibration part
    CAL <- CalibrationLasso(d,X,c=c1,maxIterPen=5e1,PostLasso=T,trace=F,maxIter=1e6)
    # W <- (1-d) * exp(X%*%CAL$betaPL)/sum(d)
    # BC <- t(d/sum(d) - W)%*%X
    
    ### 3. Computes the orthogonality parameter, using method WLS Lasso
    ORT_WLS_L <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="WLSLasso",
                                  c=c2, nopenset=c(1), RescaleY=F,
                                  maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,PostLasso=F,trace=F)
    
    #ORT_WLS_PL <- OrthogonalityReg(y,d,X,CAL$betaPL,method="WLSLasso",
    #                               c=c2, nopenset=c(1), RescaleY=F,
    #                               maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,PostLasso=T,trace=F)
    
    ### 4. Logit Lasso estimate
    LOGIT <- LogitLasso(d,X,c=.6,
                        maxIterPen=5e1,PostLasso=T,trace=F)
    
    ### 4 bis. Farrell (2015)
    FARRELL <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="LinearOutcome",
                                c=2, nopenset=c(1), RescaleY=F,
                                maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,PostLasso=T,trace=T)
    
    ### 5. BCH (2014) Estimate
    #BCH <- BCHDoubleSelec(y,d,X,cd=.95,cy=2,
    #                      nopenset=c(1),RescaleY=F,
    #                      maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,trace=F)
    
    ### 6. Third step: ATT estimation
    Results[r,] <- c(ImmunizedATT(y,d,X,CAL$betaLasso, Immunity=F)$theta,
                     #ImmunizedATT(y,d,X,CAL$betaPL, Immunity=F)$theta,
                     ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$theta,
                     #ImmunizedATT(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL, Immunity=T)$theta,
                     ImmunizedATT(y,d,X,LOGIT$betaLasso, Immunity=F)$theta,
                     #ImmunizedATT(y,d,X,LOGIT$betaPL, Immunity=F)$theta,
                     #BCH$theta,
                     ImmunizedATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso, Immunity=T)$theta,
                     ImmunizedATT(y,d,X,LOGIT$betaPL, FARRELL$muPL, Immunity=T)$theta)
    
    AsySD[r,] <- c(ImmunizedATT(y,d,X,CAL$betaLasso, Immunity=F)$sigma,
                   #ImmunizedATT(y,d,X,CAL$betaPL, Immunity=F)$sigma,
                   ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$sigma,
                   #ImmunizedATT(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL, Immunity=T)$sigma,
                   ImmunizedATT(y,d,X,LOGIT$betaLasso, Immunity=F)$sigma,
                   #ImmunizedATT(y,d,X,LOGIT$betaPL, Immunity=F)$sigma,
                   #BCH$sigma,
                   ImmunizedATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso, Immunity=T)$sigma,
                   ImmunizedATT(y,d,X,LOGIT$betaPL, FARRELL$muPL, Immunity=T)$sigma)
    
    Convergence[r,] <- c(CAL$convergence,
                         ORT_WLS_L$convergence
				#,
                        #ORT_WLS_PL$convergence
				)
    
    setTxtProgressBar(pb, r/R)
  }
  
  close(pb)
  tps <- Sys.time()-t_start
  print(tps)
  
  ## STEP B. POST-SIMULATION TREATMENT
  
  # Post-simulation treatment
  # Discard all draws which did not converge
  valid <- Convergence[,1] == 0 & Convergence[,2]==0  #& Convergence[,3] == 0
  Res <- Results[valid,]
  Res <- Res[!is.na(Res[,1]),]
  R <- nrow(Res)
  
  # Draw the charts
  id <- c(mapply(function(x) rep(x,R),1:5))
  val <- c(Res)
  data_res <- data.frame(val = val, model = id)
  
  M <- max(abs(quantile(Res,.01)),abs(quantile(Res,.99)))
  lb <- -1.1*M
  ub <- 1.1*M
  
  
  sdBCH <- mean(AsySD[,2])
  
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
  
  
  pdf(paste("plots_simu/PlugInvsImmunized_",N,"_",P,".pdf", sep=""), width=12, height=4)
  grid.arrange(get.plot(data_res,1,"Naive Plug-In, Lasso", sdBCH), 
              get.plot(data_res,2,"Immunized, Lasso", sdBCH), ncol=2)
  dev.off()
  
  pdf(paste("plots_simu/LogitFarrell_",N,"_",P,".pdf", sep=""), width=12, height=12)
  grid.arrange(get.plot(data_res,3,"IPW Logit, Lasso", sdBCH), 
               #get.plot(data_res,6,"IPW Logit, Post-Lasso", sdBCH),
               get.plot(data_res,4,"Farrell, Lasso", sdBCH), 
               get.plot(data_res,5,"Farrell, Post-Lasso", sdBCH),
               ncol=2,nrow=2)
  dev.off()
  
  ### Compute bias and RMSE
  StatDisplay <- data.frame()
  StatDisplay[1:5,"bias"] <- apply(Res,2,mean)
  StatDisplay[1:5,"RMSE"]  <- sqrt(apply(Res^2,2,mean))
  StatDisplay[1:5,"AsySD"]  <- apply(AsySD,2,mean)
  StatDisplay[1:5,"ShapiroTest"]  <- apply(Res,2, function(x) shapiro.test(x)$p.value)
  row.names(StatDisplay) <- c("NPILasso",  #"NPIPL",
                              "ImmunizedLasso",#"ImmunizedPL",
                              "LogitLasso", #"LogitPL", 
                              #"BCH 2014", 
                              "Farrell, Lasso", "Farrell, Post-Lasso")
  #print(StatDisplay)
  
  
  return(list(Results, AsySD, Convergence, StatDisplay))
  
}


#########
#########


N50P50 <- Simu_heter(N=50,P=50)
N50P100 <- Simu_heter(N=50,P=100)

N100P50 <- Simu_heter(N=100,P=50)
N100P100 <- Simu_heter(N=100,P=100)
N100P200 <- Simu_heter(N=100,P=200)

N200P50 <- Simu_heter(N=200,P=50)
N200P100 <- Simu_heter(N=200,P=100)
N200P200 <- Simu_heter(N=200,P=200)
N200P500 <- Simu_heter(N=200,P=500)

N500P50 <- Simu_heter(N=500,P=50)
N500P100 <- Simu_heter(N=500,P=100)
N500P200 <- Simu_heter(N=500,P=200)
N500P500 <- Simu_heter(N=500,P=500)


#########
#########

save.image("W:\\Bureau\\calage_temp\\BEAST\\MC_heter")

N50P50[[4]]
N50P100[[4]]

N100P50[[4]]
N100P100[[4]]
N100P200[[4]]

N200P50[[4]]
N200P100[[4]]
N200P200[[4]]
N200P500[[4]]

N500P50[[4]]
N500P100[[4]]
N500P200[[4]]
N500P500[[4]]




P50 <- rbind(N50P50[[4]],N100P50[[4]],N200P50[[4]],N500P50[[4]])
P100 <- rbind(N50P100[[4]],N100P100[[4]],N200P100[[4]],N500P100[[4]])
P200 <- rbind(N200P200[[4]],N500P200[[4]])
P500 <- rbind(N200P500[[4]],N500P500[[4]])

res <- data.frame()

res[1:20,1] <- P50[,2]
res[1:20,2] <- P50[,1]
res[1:20,3] <- P50[,4]

res[1:20,4] <- P100[,2]
res[1:20,5] <- P100[,1]
res[1:20,6] <- P100[,4]

res[11:20,7] <- P200[,2]
res[11:20,8] <- P200[,1]
res[11:20,9] <- P200[,4]

res[11:20,10] <- P500[,2]
res[11:20,11] <- P500[,1]
res[11:20,12] <- P500[,4]

res <- round(res,digits=3)
row.names(res) <- c("Naive Plug-in","BEAST","Inv. prop. weighting","Farell","Farell PL",
	"100 Naive Plug-in","100 BEAST","100 Inv. prop. weighting","100 Farell","100 Farell PL",
	"200 Naive Plug-in","200 BEAST","200 Inv. prop. weighting","200 Farell","200 Farell PL",
	"500 Naive Plug-in","500 BEAST","500 Inv. prop. weighting","500 Farell","500 Farell PL")

names(res) <- c("RMSE","Bias","Shapiro","100 RMSE","100 Bias","100 Shapiro",
	"200 RMSE","200 Bias","200 Shapiro","500 RMSE","500 Bias","500 Shapiro")

write.table(res, file = "table_simu_heter.txt", append = FALSE, quote = FALSE, sep = " & ",
            eol = paste(" \\\\ \n"), na = "--", dec = ".", row.names = T,
            col.names = T)
