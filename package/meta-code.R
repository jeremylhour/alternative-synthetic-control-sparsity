### BEAST package meta-code
### Jeremy L'Hour
### 23 fevrier 2016

library("devtools")
devtools::install_github("klutometis/roxygen")
library(roxygen2)

setwd("R:/Simulations/R_Code/package")
create("BEAST")

setwd("./BEAST")
document()

setwd("..")
install("BEAST")