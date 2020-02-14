### Parametric Alternative to SYnthetic Control: Monte Carlo Simulation
### reports RMSE, bias and coverage rate
### also parallelized
### Jeremy L'Hour
### 14/02/2020


setwd("W:/Telechargements/BEAST-master/BEAST-master")
rm(list=ls())


##############################
##############################
### PACKAGES AND FUNCTIONS ###
##############################
##############################

### Load packages
library("MASS")
library("foreach")
library("doParallel")

### Load user-defined functions
source("functions/DataSim.R") 
source("functions/DataSim_noX.R") 
source("functions/DataSim_interaction.R")
source("functions/LassoFISTA.R") 
source("functions/CalibrationLasso.R")
source("functions/OrthogonalityReg.R")
source("functions/LogitLasso.R")
source("functions/BCHDoubleSelec.R")
source("functions/ImmunizedATT.R")

func_liste = c('DataSim','DataSim_noX','DataSim_interaction',
               'LassoFISTA','CalibrationLasso','OrthogonalityReg','LogitLasso','BCHDoubleSelec','ImmunizedATT',
               'gamma','gammagrad','prox','LeastSq','LeastSqgrad','LassoObj','Logitloss','Logitlossgrad') # list of functions for running parallel loop

### Monte Carlo Simulations -- setting up the function

Simu <- function(N,P,R=100,R2y=.8,R2d=.2,c1=.7,c2=2,Table="base"){
  print(paste('--- Simulations start : R=',R,', n=',N,', p=',P,' ---'))
  print(paste('--- DGP style :',Table))
  ## STEP A. SIMULATIONS
  cores = detectCores()
  cl = makeCluster(cores[1]) #not to overload your computer
  registerDoParallel(cl)
  
  t_start <- Sys.time()
  
  resPAR <- foreach(r = 1:R,.export=func_liste,.packages=c('lbfgs','MASS'),.combine='rbind', .multicombine=TRUE) %dopar% {
    ### 1. Generate data
    if(Table=="base"){
      data <- DataSim(n=N,p=P,Ry=R2y,Rd=R2d,TreatHeter=F)
    } else if(Table=="noX") {
      data <- DataSim_noX(n=N,p=P,Ry=R2y,Rd=R2d,TreatHeter=F)
    } else if(Table=="heterogeneous"){
      data <- DataSim(n=N,p=P,Ry=R2y,Rd=R2d,TreatHeter=T)
    } else if(Table=="interaction"){
      data <- DataSim_interaction(n=N,p=P,Ry=R2y,Rd=R2d,TreatHeter=F)
    }
    
    X=data$X; y=data$y; d=data$d
    
    ### 2. Calibration part
    CAL <- CalibrationLasso(d,X,c=c1,maxIterPen=5e1,PostLasso=T,trace=F)
    
    ### 3. Computes the orthogonality parameter, using method WLS Lasso
    ORT_WLS_L <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="WLSLasso",
                                  c=c2, nopenset=c(1), RescaleY=F,
                                  maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,PostLasso=F,trace=F)
    
    ### 4. Logit Lasso estimate
    LOGIT <- LogitLasso(d,X,c=.6,maxIterPen=5e1,PostLasso=T,trace=F)
    
    ### 4 bis. Farrell (2015)
    FARRELL <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="LinearOutcome",
                                c=2, nopenset=c(1), RescaleY=F,
                                maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,PostLasso=T,trace=T)
    
    ### 5. BCH (2014) Estimate
    BCH <- BCHDoubleSelec(y,d,X,cd=.95,cy=2,
                          nopenset=c(1),RescaleY=F,
                          maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,trace=F)
    
    ### 6. Third step: ATT estimation
    Estimate <- c(ImmunizedATT(y,d,X,CAL$betaLasso, Immunity=F)$theta,
                  ImmunizedATT(y,d,X,LOGIT$betaLasso, Immunity=F)$theta,
                  ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$theta,
                  BCH$theta,
                  ImmunizedATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso, Immunity=T)$theta,
                  ImmunizedATT(y,d,X,LOGIT$betaPL, FARRELL$muPL, Immunity=T)$theta)
    
    AsySD <- c(ImmunizedATT(y,d,X,CAL$betaLasso, Immunity=F)$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaLasso, Immunity=F)$sigma,
               ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$sigma,
               BCH$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso, Immunity=T)$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaPL, FARRELL$muPL, Immunity=T)$sigma)
    
    
    Convergence <- c(CAL$convergence,
                         ORT_WLS_L$convergence)
    
    c(Estimate,AsySD,Convergence)
  }
  
  print('--- Simulations over ! ---')
  print(paste('--- Time Elapsed : ',Sys.time()-t_start,' ---'))
  stopCluster(cl)
  
  Estimate = resPAR[,1:6]; AsySD = resPAR[,7:12]; Convergence =  resPAR[,13:14]
  
  ## STEP B. POST-SIMULATION TREATMENT
  
  # Post-simulation treatment
  # Discard all draws which did not converge
  valid = (Convergence[,1] == 0 & Convergence[,2]==0)
  Estimate = Estimate[valid,]
  AsySD = AsySD[valid,]
  
  ### Compute bias and RMSE
  StatDisplay = data.frame()
  StatDisplay[1:6,"RMSE"]  = sqrt(apply(Estimate^2,2,mean))
  StatDisplay[1:6,"bias"] = apply(Estimate,2,mean)
  
  borne_sup = Estimate + qnorm(.975) * AsySD
  borne_inf = Estimate - qnorm(.975) * AsySD
  StatDisplay[1:6,"CoverageRate"] = apply((borne_sup>0)*(borne_inf<0),2,mean) # zero if zero is not included in CI
  
  
  row.names(StatDisplay) <- c("Naive Lasso","IPW LogitLasso","Immunized Lasso",
                              "BCH 2014","Farrell Lasso","Farrell Post-Lasso")
  print(StatDisplay)
  
  return(list(Estimate = Estimate, AsySD = AsySD, StatDisplay = StatDisplay))
}


###########################
###########################
### RUNNING SIMULATIONS ###
###########################
###########################

DGP_style = 'noX' # to modify to generate each table

set.seed(30031987)

N50P50 <- Simu(N=50,P=50,Table=DGP_style)
N50P100 <- Simu(N=50,P=100,Table=DGP_style)

N100P50 <- Simu(N=100,P=50,Table=DGP_style)
N100P100 <- Simu(N=100,P=100,Table=DGP_style)
N100P200 <- Simu(N=100,P=200,Table=DGP_style)

N200P50 <- Simu(N=200,P=50,Table=DGP_style)
N200P100 <- Simu(N=200,P=100,Table=DGP_style)
N200P200 <- Simu(N=200,P=200,Table=DGP_style)
N200P500 <- Simu(N=200,P=500,Table=DGP_style)

N500P50 <- Simu(N=500,P=50,Table=DGP_style)
N500P100 <- Simu(N=500,P=100,Table=DGP_style)
N500P200 <- Simu(N=500,P=200,Table=DGP_style)
N500P500 <- Simu(N=500,P=500,Table=DGP_style)


#########################
#########################
### FORMATING RESULTS ###
#########################
#########################

save.image("W:\\Bureau\\calage_temp\\BEAST\\MC_base")


P50 <- rbind(N50P50$StatDisplay,N100P50$StatDisplay,N200P50$StatDisplay,N500P50$StatDisplay)
P100 <- rbind(N50P100$StatDisplay,N100P100$StatDisplay,N200P100$StatDisplay,N500P100$StatDisplay)
P200 <- rbind(N200P200$StatDisplay,N500P200$StatDisplay)
P500 <- rbind(N200P500$StatDisplay,N500P500$StatDisplay)

res <- data.frame()

res[1:24,1:3] <- P50
res[1:24,4:6] <- P100
res[13:24,7:9] <- P200
res[13:24,10:12] <- P500

res <- round(res,digits=3)

estim_names <- c("Naive Plug-in","IPW Logit-Lasso","Immunized","BCH","Farell","Farell PL")
row.names(res) <- c(paste('50',estim_names),
                    paste('100',estim_names),
                    paste('200',estim_names),
                    paste('500',estim_names))

names(res) <- c("RMSE","Bias","Coverage Rate","RMSE","Bias","Coverage Rate",
                "RMSE","Bias","Coverage Rate","RMSE","Bias","Coverage Rate")

write.table(res, file = paste("Table_",DGP_style,".txt",sep=''), append = FALSE, quote = FALSE, sep = " & ",
            eol = paste(" \\\\ \n"), na = "--", dec = ".", row.names = T,
            col.names = T)