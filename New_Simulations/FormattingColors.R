#######################################
#######################################
### Mise en couleur des résultats #####
#######################################
#######################################
library("data.table") # loads the smart "fast" fread function
rm(list=ls())

#import data into R
DGP_style="noX"
data <- fread("/Users/jeremylhour/Downloads/output/Table_noX.txt", skip=1, sep="&", data.table=F) # skips in compatible lines
data <- data[,-1]
names <- read.table("/Users/jeremylhour/Downloads/output/Table_noX.txt", skip=0, nrow=1, sep="&", stringsAsFactors=F)

names = gsub("\\\\", "", as.character(names))
names = gsub(" ", "", as.character(names))
colnames(data) = names

data[,ncol(data)] = gsub("\\\\", "", as.character(data[,ncol(data)]))
data[,ncol(data)] = gsub(" ", "", as.character(data[,ncol(data)]))
data[,ncol(data)] = as.numeric(data[,ncol(data)]) #change type as desired

data[data=="--"] <- as.numeric(NA)
for(k in 1:ncol(data)){
  data[,k] = as.numeric(data[,k])
}
res = data

####

nb_e = 7

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

estim_names <- c('Naive Plug-in','IPW Logit-Lasso',
                 'Immunized','Immunized Unweighted Reg',
                 'BCH','Farell','Farell PL')

row.names(res_colored) <- c(paste('n=50',estim_names),
                            paste('n=100',estim_names),
                            paste('n=200',estim_names),
                            paste('n=500',estim_names))

write.table(res_colored, file = paste("/Users/jeremylhour/Downloads/output/ColorTable_",DGP_style,".txt",sep=''), append = FALSE, quote = FALSE, sep = " & ",
            eol = paste(" \\\\ \n"), na = "--", dec = ".", row.names = T,
            col.names = T)