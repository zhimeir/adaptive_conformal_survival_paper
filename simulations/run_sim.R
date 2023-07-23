########################################
## Process the input argument
########################################
args <- commandArgs(trailingOnly = TRUE)
seed<- as.integer(args[1])
if(is.na(seed)){seed <- 1}


########################################
## load libraries
########################################
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(quantreg))
suppressPackageStartupMessages(library(GauPro))
suppressPackageStartupMessages(library(gbm))
suppressPackageStartupMessages(library(grf))
## suppressPackageStartupMessages(library(caret))

########################################
### source code
########################################
source("./source_code.R")
source("./model_script.R")
source("./simu.R")


########################################
### run simulations
########################################
## configurations
setting_list = c("ld_setting1",
                 "ld_setting2",
                 "ld_setting3",
                 "ld_setting4",
                 "hd_homosc",
                 "hd_heterosc")

## parameters
alpha <- .1    # target level 1-alpha
n <- 1000
n_test <- 5000
n_train <- n
n_calib <- n
xmin <- 0 
xmax <- 4
beta <- 20 / sqrt(n)
exp_rate <- .1

for(setting in setting_list){
  if(setting %in% c("hd_homosc","hd_heterosc")){
    p <- 10
  }else{
    p <- 1
  }
  simures <- simu(seed + 1234, setting, n, p, 
                   n_train, n_calib, n_test, 
                   beta, xmin, xmax, 
                   exp_rate, alpha, "cox") 
   save_dir <- sprintf("../results/%s_seed_%d.csv", setting, seed)
   write.csv(simures, save_dir)
}
