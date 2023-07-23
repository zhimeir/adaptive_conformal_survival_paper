########################################
## load libraries
########################################
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(survival))
suppressPackageStartupMessages(library(quantreg))
suppressPackageStartupMessages(library(conTree))
suppressPackageStartupMessages(library(GauPro))
suppressPackageStartupMessages(library(gbm))
suppressPackageStartupMessages(library(grf))
suppressPackageStartupMessages(library(caret))
library(foreach)
library(doParallel)
library(ggplot2)
library(ggExtra)
source("~/Documents/GitHub/adaptive_cfsurv-main/workflow/simulation/pings/pings_run.R")
datapath = "~/Documents/GitHub/adaptive_cfsurv-main/workflow/simulation/pings/"

m0 = 14
L = 21
K = 16 # number of dataset splitting


data <- read.csv(paste0(datapath,"data_", m0,".csv"))
data = data[,-1]


N0 = nrow(data)
m = 0
seed = 1234 + 7
set.seed(seed)
samp = NULL
n_test = 500
n = 1000+n_test
n_train = (n-n_test)/2
tr = 1:n_train
dat_cp = setdiff(1:N0,tr)
for(j in 1:K){
  samp = rbind(samp,c(tr,sample(dat_cp,(n+n_test)/2)))
}
T1 = data$censored_T
trueT = T1
c_list = NULL



detectCores()
cl<- makeCluster(detectCores(), outfile="")
registerDoParallel(cl)
alpha.seq = c(0.1)
nam = c("adaptive LPB (T)","adaptive LPB (CT)","Fix-LPB","Naive-CQR","Cox","RF")
num_method = length(nam)
mod_list = c("cox")
for(alpha in alpha.seq){
  print(alpha)
  simures = foreach(i = 1:K, .combine="rbind"
  ) %dopar% real.simu(data[samp[i,],],trueT[samp[i,]],alpha,n_test,c_list)
  write.csv(simures,file = paste0("~/Documents/GitHub/adaptive_cfsurv-main/result/ping/hist_",
                                  m0,"_",m,"_",seed,"_",K,"_",n,"_",n_test,"_",alpha*100,"_",opt,".csv"))
}

stopCluster(cl) 
