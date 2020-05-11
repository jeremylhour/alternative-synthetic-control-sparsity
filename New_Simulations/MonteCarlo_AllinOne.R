### Parametric Alternative to Synthetic Control: Monte Carlo Simulation
### reports RMSE, bias and coverage rate, also parallelized
### Jeremy L'Hour
### 14/02/2020
### Last edited: 13/03/2020


setwd("W:/1A_These/A. Research/beast_git/BEAST")
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
library('lbfgs')

### Load user-defined functions
source("functions/DataSim.R") 
source("functions/New_DataSim.R") 
source("functions/DataSim_New2.R") 
source("functions/DataSim_noX.R") 
source("functions/DataSim_interaction.R")
source("functions/LassoFISTA.R") 
source("functions/CalibrationLasso.R")
source("functions/OrthogonalityReg.R")
source("functions/LogitLasso.R")
source("functions/BCHDoubleSelec.R")
source("functions/ImmunizedATT.R")

func_liste = c('DataSim','DataSim_noX','DataSim_interaction','New_DataSim','sigmoid', 'DataSim_New2',
               'LassoFISTA','CalibrationLasso','OrthogonalityReg','LogitLasso','BCHDoubleSelec','ImmunizedATT',
               'gamma','gammagrad','prox','LeastSq','LeastSqgrad','LassoObj','Logitloss','Logitlossgrad') # list of functions for running parallel loop

### Monte Carlo Simulations -- setting up the function

Simu <- function(N,P,R=10000,R2y=.8,R2d=.3,Table="base"){
  print(paste('--- Simulations start : R =',R,', n =',N,', p =',P,' ---'))
  print(paste('--- DGP style :',Table,' ---'))
  ## STEP A. SIMULATIONS
  cores = detectCores()
  cl = makeCluster(20) #not to overload your computer
  registerDoParallel(cl)
  
  t_start <- Sys.time()
  
  resPAR <- foreach(r = 1:R,.export=func_liste,.packages=c('lbfgs','MASS'),.combine='rbind', .multicombine=TRUE, .errorhandling = 'remove') %dopar% {
    ### 1. Generate data
    ATT = 0
    if(Table=="base"){
      data <- DataSim(n=N,p=P,Ry=R2y,Rd=R2d,TreatHeter=F)
    } else if(Table=="noX") {
      data <- DataSim_noX(n=N,p=P,Ry=R2y,Rd=R2d,TreatHeter=F)
    } else if(Table=="heterogeneous"){
      data <- DataSim(n=N,p=P,Ry=R2y,Rd=R2d,TreatHeter=T)
    } else if(Table=="interaction"){
      data <- DataSim_interaction(n=N,p=P,Ry=R2y,Rd=R2d,TreatHeter=F)
    } else if(Table=="newdgp"){
      data <- New_DataSim(n=N,p=P,Ry=R2y,Rd=R2d)
      ATT <- data$ATT
    } else if(Table=="newdgp2"){
      data <- DataSim_New2(n=N,p=P,Ry=R2y,Rd=R2d)
      ATT <- data$ATT
    }
    
    X=data$X; y=data$y; d=data$d
    
    ### 2. Calibration part
    CAL <- CalibrationLasso(d,X,c=1.001,maxIterPen=5e1,PostLasso=T,trace=F)
    
    ### 3. Computes the orthogonality parameter, using method WLS Lasso
    ORT_WLS_L <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="WLSLasso",
                                           c=1.001, nopenset=c(1), RescaleY=F,
                                           maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,PostLasso=F,trace=F)
    
    ### 3bis. Computes the orthogonality parameter, using method WLS Post-Lasso
    ORT_WLS_PL <- OrthogonalityReg(y,d,X,CAL$betaPL,method="WLSLasso",
                                  c=1.001, nopenset=c(1), RescaleY=F,
                                  maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,PostLasso=T,trace=F)
    
    ### 4. Logit Lasso estimate
    LOGIT <- LogitLasso(d,X,c=1.001,maxIterPen=5e1,PostLasso=T,trace=F)
    
    ### 4 bis. Farrell (2015)
    FARRELL <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="LinearOutcome",
                                c=1.001, nopenset=c(1), RescaleY=F,
                                maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,PostLasso=T,trace=T)
    
    ### 5. BCH (2014) Estimate
    BCH <- BCHDoubleSelec(y,d,X,cd=1.001,cy=1.001,
                          nopenset=c(1),RescaleY=F,
                          maxIterPen=1e4,maxIterLasso=1e4,tolLasso=1e-6,trace=F)
    
    ### 6. Naive Oracle (for New_DGP)
    ORACLE <- CalibrationLasso(d,X[,c(1:11,(P-8):(P+1))],c=0,maxIterPen=5e1,PostLasso=F,trace=F)
    
    ### 7. Third step: ATT estimation
    Estimate <- c(ImmunizedATT(y,d,X,CAL$betaLasso, Immunity=F)$theta,
                  ImmunizedATT(y,d,X,LOGIT$betaLasso, Immunity=F)$theta,
                  ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$theta,
                  ImmunizedATT(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL, Immunity=T)$theta,
                  BCH$theta,
                  ImmunizedATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso, Immunity=T)$theta,
                  ImmunizedATT(y,d,X,LOGIT$betaPL, FARRELL$muPL, Immunity=T)$theta,
                  ImmunizedATT(y,d,X[,c(1:11,(P-8):(P+1))],ORACLE$betaLasso, Immunity=F)$theta)
    
    AsySD <- c(ImmunizedATT(y,d,X,CAL$betaLasso, Immunity=F)$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaLasso, Immunity=F)$sigma,
               ImmunizedATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso, Immunity=T)$sigma,
               ImmunizedATT(y,d,X,CAL$betaPL,ORT_WLS_PL$muPL, Immunity=T)$sigma,
               BCH$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaLasso, FARRELL$muLasso, Immunity=T)$sigma,
               ImmunizedATT(y,d,X,LOGIT$betaPL, FARRELL$muPL, Immunity=T)$sigma,
               ImmunizedATT(y,d,X[,c(1:11,(P-8):(P+1))],ORACLE$betaLasso, Immunity=F)$sigma)
    
    
    Convergence <- c(CAL$convergence,
                         ORT_WLS_L$convergence)
    
    c(Estimate,AsySD,Convergence,ATT)
  }
  
  print('--- Simulations over ! ---')
  print(Sys.time()-t_start)
  stopCluster(cl)
  
  nb_e = 8 # nb of estimators
  
  Estimate = resPAR[,1:nb_e]; AsySD = resPAR[,(nb_e+1):(2*nb_e)]; Convergence =  resPAR[,(2*nb_e+1):((2*nb_e+2))]
  ATT = mean(resPAR[,ncol(resPAR)])
  #ATT = resPAR[,ncol(resPAR)]
  
  ## STEP B. POST-SIMULATION TREATMENT
  
  # Post-simulation treatment
  # Discard all draws which did not converge
  valid = (Convergence[,1] == 0 & Convergence[,2]==0)
  Estimate = Estimate[valid,]
  AsySD = AsySD[valid,]
  
  print(paste('--- Convergence pct : ',round(sum(valid)/R,digits=3),' ---'))
  
  
  ### Compute bias, RMSE and Coverage Rate
  StatDisplay = data.frame()
  StatDisplay[1:nb_e,"RMSE"]  = sqrt(apply((Estimate-ATT)^2,2,mean,na.rm=T))
  StatDisplay[1:nb_e,"bias"] = apply((Estimate-ATT),2,mean,na.rm=T)
  
  borne_sup = Estimate + qnorm(.975) * AsySD - ATT
  borne_inf = Estimate - qnorm(.975) * AsySD - ATT
  StatDisplay[1:nb_e,"CoverageRate"] = apply((borne_sup>0)*(borne_inf<0),2,mean,na.rm=T) # zero if zero is not included in CI
  
  
  row.names(StatDisplay) <- c("Naive -- Balancing Lasso","Naive -- IPW Logit Lasso","Immunized -- Lasso",
                              "Immunized -- Post-Lasso","BCH 2014","Farrell Lasso","Farrell Post-Lasso",
                              "Oracle -- Balancing")
  print(round(StatDisplay,digits=3))
  
  return(list(Estimate = Estimate, AsySD = AsySD, StatDisplay = StatDisplay))
}


###########################
###########################
### RUNNING SIMULATIONS ###
###########################
###########################

DGP_style = "newdgp" # modify here to generate each table

set.seed(9081993)

N500P50 <- Simu(N=500,P=50,Table=DGP_style)
N500P100 <- Simu(N=500,P=100,Table=DGP_style)
N500P200 <- Simu(N=500,P=200,Table=DGP_style)
N500P500 <- Simu(N=500,P=500,Table=DGP_style)

N200P50 <- Simu(N=200,P=50,Table=DGP_style)
N200P100 <- Simu(N=200,P=100,Table=DGP_style)
N200P200 <- Simu(N=200,P=200,Table=DGP_style)
N200P500 <- Simu(N=200,P=500,Table=DGP_style)

N100P50 <- Simu(N=100,P=50,Table=DGP_style)
N100P100 <- Simu(N=100,P=100,Table=DGP_style)
N100P200 <- Simu(N=100,P=200,Table=DGP_style)

N50P50 <- Simu(N=50,P=50,Table=DGP_style)
N50P100 <- Simu(N=50,P=100,Table=DGP_style)


#save.image("//ulysse/users/JL.HOUR/1A_These/sim_output")

#########################
#########################
### FORMATING RESULTS ###
#########################
#########################

estim_names <- c("Naive Plug-In -- Balancing Lasso","Naive Plug-In -- Logit Lasso","Immunized -- Lasso",
                 "Immunized -- Post-Lasso","BCH 2014","Farrell Lasso","Farrell Post-Lasso",
                 "Oracle -- Balancing")
nb_e = length(estim_names)

res <- data.frame()

res[1:(4*nb_e),1:3] <- rbind(N50P50$StatDisplay,N100P50$StatDisplay,N200P50$StatDisplay,N500P50$StatDisplay) # p = 50
res[1:(4*nb_e),4:6] <- rbind(N50P100$StatDisplay,N100P100$StatDisplay,N200P100$StatDisplay,N500P100$StatDisplay) # p = 100
res[(2*nb_e+1):(4*nb_e),7:9] <- rbind(N200P200$StatDisplay,N500P200$StatDisplay) # p = 200
res[(2*nb_e+1):(4*nb_e),10:12] <- rbind(N200P500$StatDisplay,N500P500$StatDisplay) # p = 500

res <- round(res,digits=3)

row.names(res) <- c(paste('n=50',estim_names),
                    paste('n=100',estim_names),
                    paste('n=200',estim_names),
                    paste('n=500',estim_names))

names(res) <- rep(c("RMSE","Bias","Cov. Rate"),4)


write.table(res, file = paste("New_Simulations/Table_",DGP_style,".txt",sep=''), append = FALSE, quote = FALSE, sep = " & ",
            eol = paste(" \\\\ \n"), na = "--", dec = ".", row.names = T,
            col.names = T)

####################################
####################################
### FUNCTIONS TO COLOR RESULTS #####
####################################
####################################

mm_rescale <- function(x,absolute=F){
  if(absolute){
    y = abs(x)
  } else {
    y = x
  }
  
  if(max(y)-min(y) == 0){
    y_rescaled = round(y-min(y),digits=3)
  } else {
    y_rescaled = round((y-min(y))/(max(y)-min(y)),digits=3)
  }
  
  return(y_rescaled)
}

col_latex <- function(rouge,vert,bleu){
  paste("\\cellcolor[rgb]{",rouge,",",vert,",",bleu,"}",sep='') 
}

return_col <- function(x,shading='bleu',absolute=F){
  y = mm_rescale(x,absolute)
  if(shading=='bleu'){
    col_text = mapply(function(x) col_latex(.31+.69*x,.47+x/2,bleu = .753),y)
  } else if (shading=='rouge'){
    col_text = mapply(function(x) col_latex(rouge=1,.95-x/2,.55-x/2),y)
  }
  col_text = paste(col_text,round(x,digits=3))
  return(col_text)
}

#####################
#####################
### COLOR RESULTS ###
#####################
#####################

res_colored = res

for(bloc in 1:4){
  for(k in 1:ncol(res)){
    if(sum(is.na(res_colored[((bloc-1)*nb_e+1):(bloc*nb_e),k]))>0){
      next 
    }
    
    if(k %% 3==1){
      res_colored[((bloc-1)*nb_e+1):(bloc*nb_e),k] = return_col(res[((bloc-1)*nb_e+1):(bloc*nb_e),k],shading="bleu")
    } else if(k %% 3==2){
      res_colored[((bloc-1)*nb_e+1):(bloc*nb_e),k] = return_col(res[((bloc-1)*nb_e+1):(bloc*nb_e),k],shading="bleu",absolute=T)
    } else if(k %% 3==0){
      res_colored[((bloc-1)*nb_e+1):(bloc*nb_e),k] = return_col(res[((bloc-1)*nb_e+1):(bloc*nb_e),k],shading="rouge")
    }
  }
}

write.table(res_colored, file = paste("New_Simulations/ColorTable_",DGP_style,".txt",sep=''), append = FALSE, quote = FALSE, sep = " & ",
            eol = paste(" \\\\ \n"), na = "--", dec = ".", row.names = T,
            col.names = T)