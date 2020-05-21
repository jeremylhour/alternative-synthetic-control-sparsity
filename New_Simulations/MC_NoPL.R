### Parametric Alternative to Synthetic Control: Monte Carlo Simulation
### reports RMSE, bias and coverage rate, also parallelized
### Jeremy L'Hour
### 14/02/2020
### Last edited: 14/05/2020

### NO POST-LASSO


setwd("W:/1A_These/A. Research/beast_git/BEAST")
#setwd("/Users/jeremylhour/Documents/BEAST")
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
source("functions/Compute_ATT.R")

func_liste = c('DataSim','DataSim_noX','DataSim_interaction','New_DataSim','sigmoid', 'DataSim_New2',
               'LassoFISTA','CalibrationLasso','OrthogonalityReg','LogitLasso','BCHDoubleSelec','Compute_ATT',
               'gamma','gammagrad','prox','LeastSq','LeastSqgrad','LassoObj','Logitloss','Logitlossgrad','sigmoid') # list of functions for running parallel loop

### Monte Carlo Simulations -- setting up the function

Simu <- function(N,P,R=100,R2y=.8,R2d=.3,Table="base"){
  print(paste('--- Simulations start : R =',R,', n =',N,', p =',P,' ---'))
  print(paste('--- DGP style :',Table,' ---'))
  ## STEP A. SIMULATIONS
  cores = detectCores()
  cl = makeCluster(7,setup_timeout=.5) #not to overload your computer
  registerDoParallel(cl)
  
  t_start <- Sys.time()
  
  resPAR <- foreach(r = 1:R,.export=func_liste,.packages=c('lbfgs','MASS'),.combine='rbind', .multicombine=TRUE, .errorhandling = 'remove') %dopar% {
    set.seed(r);
    ### 0. Generate data
    if(Table=="newdgp"){
      data <- New_DataSim(n=N,p=P,Ry=R2y,Rd=R2d)
    } else if(Table=="newdgp2"){
      data <- DataSim_New2(n=N,p=P,Ry=R2y,Rd=R2d)
    }
    ATT <- data$ATT
    X=data$X; y=data$y; d=data$d
    
    max_iter_penalty = 300 # can be important in term of computation time
    
    ### 1. Compute Balancing parameter
    CAL <- CalibrationLasso(d,X,c=1.001,maxIterPen=max_iter_penalty,PostLasso=F,trace=F)
    
    ### 2. Computes the orthogonality parameter, using method WLS Lasso
    ORT_WLS_L <- OrthogonalityReg(y,d,X,CAL$betaLasso,method="WLSLasso",
                                  c=1.001, nopenset=c(1), RescaleY=F,
                                  maxIterPen=max_iter_penalty,maxIterLasso=1e4,tolLasso=1e-6,PostLasso=F,trace=F)
    
    ### 3. Logit Lasso
    LOGIT <- LogitLasso(d,X,c=1.001,maxIterPen=max_iter_penalty,PostLasso=F,trace=F)
    
    ### 4bis. Orthogonality param for Farrell (2015)
    FARRELL <- OrthogonalityReg(y,d,X,LOGIT$betaLasso,method="LinearOutcome",
                                c=1.001, nopenset=c(1), RescaleY=F,
                                maxIterPen=max_iter_penalty,maxIterLasso=1e4,tolLasso=1e-6,PostLasso=F,trace=F)
    
    ### 5. BCH (2014) double post-selection
    BCH <- BCHDoubleSelec(y,d,X,cd=1.001,cy=1.001,
                          nopenset=c(1),RescaleY=F,
                          maxIterPen=max_iter_penalty,maxIterLasso=1e4,tolLasso=1e-6,trace=F)
    
    ### 6. Oracle
    ORACLE <- CalibrationLasso(d,X[,c(1:11)],c=0,maxIterPen=max_iter_penalty,PostLasso=F,trace=F)
    
    ### 7. Third step: ATT estimation
    Estimate <- c(Compute_ATT(y,d,X,CAL$betaLasso)$theta,
                  Compute_ATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso)$theta,
                  BCH$theta,
                  Compute_ATT(y,d,X,LOGIT$betaLasso,FARRELL$muLasso)$theta,
                  Compute_ATT(y,d,X[,c(1:11)],ORACLE$betaLasso)$theta)
    
    AsySD <- c(Compute_ATT(y,d,X,CAL$betaLasso)$sigma,
               Compute_ATT(y,d,X,CAL$betaLasso,ORT_WLS_L$muLasso)$sigma,
               BCH$sigma,
               Compute_ATT(y,d,X,LOGIT$betaLasso,FARRELL$muLasso)$sigma,
               Compute_ATT(y,d,X[,c(1:11)],ORACLE$betaLasso)$sigma)
    
    Convergence <- c(CAL$convergence,
                     ORT_WLS_L$convergence)
    
    c(Estimate,AsySD,Convergence,ATT)
  }
  
  print('--- Simulations over ! ---')
  print(Sys.time()-t_start)
  stopCluster(cl)
  
  nb_e = 5 # nb of estimators
  
  Estimate = resPAR[,1:nb_e]; AsySD = resPAR[,(nb_e+1):(2*nb_e)]; Convergence = resPAR[,(2*nb_e+1):((2*nb_e+2))]
  ATT = mean(resPAR[,ncol(resPAR)])
  
  ## STEP B. POST-SIMULATION TREATMENT
  
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
  
  
  row.names(StatDisplay) <- c("Naive -- Balancing Lasso","Immunized -- Lasso",
                              "BCH (2014)","Farrell (2015) -- Lasso",
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

# P = 50
N500P50 <- Simu(N=500,P=50,Table=DGP_style)
N1000P50 <- Simu(N=1000,P=50,Table=DGP_style)

# P = 200
N500P200 <- Simu(N=500,P=200,Table=DGP_style)
N1000P200 <- Simu(N=1000,P=200,Table=DGP_style)

# P = 500
N500P500 <- Simu(N=500,P=500,Table=DGP_style)
N1000P500 <- Simu(N=1000,P=500,Table=DGP_style)

# P = 1000
N500P1000 <- Simu(N=500,P=1000,Table=DGP_style)
N1000P1000 <- Simu(N=1000,P=1000,Table=DGP_style)


#save.image("//ulysse/users/JL.HOUR/1A_These/sim_output")

#########################
#########################
### FORMATING RESULTS ###
#########################
#########################

estim_names <- c("Naive -- Balancing Lasso","Immunized -- Lasso",
                 "BCH (2014)","Farrell (2015) -- Lasso",
                 "Oracle -- Balancing")
nb_e = length(estim_names)

res <- data.frame()

res[1:(2*nb_e),1:3] <- rbind(N500P50$StatDisplay,
                             N1000P50$StatDisplay) # p = 50
res[1:(2*nb_e),4:6] <- rbind(N500P200$StatDisplay,
                             N1000P200$StatDisplay) # p = 200
res[1:(2*nb_e),7:9] <- rbind(N500P500$StatDisplay,
                             N1000P500$StatDisplay) # p = 500
res[1:(2*nb_e),10:12] <- rbind(N500P1000$StatDisplay,
                               N1000P1000$StatDisplay) # p = 1000

res <- round(res,digits=3)

row.names(res) <- c(paste('n=500',estim_names),
                    paste('n=1000',estim_names))

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

for(bloc in 1:5){
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
