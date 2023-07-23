sigma_x <- function(x) (5 + x)/10
z = rbinom(1,1,.5)
# gen_t <- function(x) exp(2 + z * x + (1-z) * abs(xmax-x)  
#                          +  (z * sigma_x(x) + (1-z) * sigma_x(xmax-x)) * rnorm(length(x)))
# gen_c <- function(x) exp(2 + 2*x +  max(sigma_x(xmax-x),sigma_x(x)) * rnorm(length(x)))
gen_t <- function(x) exp(2 + sqrt(x) + sin(5*x) +  sigma_x(x) * rnorm(length(x)))
gen_c <- function(x) exp(2 + sqrt(x) +  sigma_x(x) * rnorm(length(x)))
# sigma_x <- function(x) (5 + x)/5
# gen_t <- function(x) exp(2 + (20/sqrt(3000)) * sqrt(abs(x)) +  sigma_x(x) * rnorm(length(x))) 
# gen_c <- function(x) rexp(rate = 0.4, n = length(x)) 

model_generating_fun <- function(n_train, n_calib, n_test, setting, beta, xnames, xmin, xmax, exp_rate){
  if(setting == "ld_heterosc20"){
    p = 1
    # sigma_x <- function(x) (5 + x)/2
    # sigma_x <- function(x) (5 + x)/20
    # ########################################
    # ## Data generating models
    # ########################################
    # gen_t <- function(x) exp(2 + beta * sqrt(abs(x)) +  sigma_x(x) * rnorm(length(x)))
    # gen_c <- function(x) exp(2 + beta * sqrt(abs(x-xmax)) +  sigma_x(x-xmax) * rnorm(length(x)))
    # sigma_x <- function(x) (5 + x)/5
    # gen_t <- function(x) exp(2 + beta * x^2 +  sigma_x(x) * rnorm(length(x)))
    # gen_c <- function(x) exp(2 + beta * (xmax-x)^2 +  sigma_x(x-xmax) * rnorm(length(x)))
    z = rbinom(1,1,.5)
    ########################################
    ## Generate training data
    ########################################
    # set.seed(24601)
    X <- runif(n_train, xmin, xmax)
    T <- gen_t(X)
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data_fit <- data.frame(X1 = X, C = C, censored_T = censored_T, event = event)
    
    ########################################
    ## Generate the calibration data and the test data
    ########################################
    # set.seed(seed)
    X <- runif(n_calib + n_test, xmin, xmax)     
    T <- gen_t(X) 
    C <- gen_c(X)
    event <- (T < C)
    censored_T <- pmin(T, C)
    data <- data.frame(X1 = X, C = C, event = event, censored_T = censored_T)
    data_calib <- data[1 : n_calib, ]
    data_test <- data[(n_calib + 1) : (n_calib + n_test),]
    data <- rbind(data_fit, data_calib)
    T_test = T[(n_calib + 1) : (n_calib + n_test)]
  }
  ## Collect results 
  obj <- list(data_fit = data_fit, 
              data_calib = data_calib, 
              data_test = data_test, 
              data = data, 
              T_test = T_test)
  
  return(obj)
  
}


run.simu = function(seed,alpha=.05,n){
  ########################################
  ## load libraries
  ########################################
  suppressPackageStartupMessages(library(tidyverse))
  suppressPackageStartupMessages(library(survival))
  suppressPackageStartupMessages(library(quantreg))
  suppressPackageStartupMessages(library(conTree))
  suppressPackageStartupMessages(library(GauPro))
  suppressPackageStartupMessages(library(gbm))
  suppressPackageStartupMessages(library(foreach))
  suppressPackageStartupMessages(library(doParallel))
  suppressPackageStartupMessages(library(grf))
  suppressPackageStartupMessages(library(caret))
  dir = "~/Documents/Github/adaptive_cfsurv-main/workflow/simulation"
  setwd(dir)
  source("source_code.R")

  ## Initialization
  set.seed(seed)
  p = 1
  xnames <- paste0("X",1:p)

  ## Generate data according to the setting
  n_test <- 100
  n_train <- n
  n_calib <- n
  n = n_train
  beta <- 20 / sqrt(n)
  xmin <- 0; xmax <- 4
  setting <- "ld_heterosc20"
  # data_obj <- simulation()
  data_obj <- model_generating_fun(n_train, n_calib, n_test, setting, beta, xnames, xmin, xmax)
  data_fit <- data_obj$data_fit
  data_calib <- data_obj$data_calib
  data_test <- data_obj$data_test
  data <- data_obj$data
  T_test <- data_obj$T_test

  ## Arguments to be passed to the subsequent functions
  x=data_test[,names(data_test) %in% xnames]
  Xtrain = data[,names(data) %in% xnames]
  C = data$C
  event = data$event
  time=data$censored_T
  type="quantile"
  dist= "weibull"
  I_fit = 1:n_train
  ftol=.1
  tol=.1
  n.tree=10
  fit = data_fit
  fit$C <- -fit$C
  mdl0 <- GauPro(X = as.matrix(fit[,names(fit) %in% xnames]),
                 Z = fit$C, D = p,
                 type = "Gauss")

  alpha_list <- alpha

  ## Fit the survival model
  newdata <- data.frame(x)
  colnames(newdata) <- xnames
  fmla <- as.formula(paste("Surv(censored_T, event) ~ ", paste(xnames, collapse= "+")))
  mdl <- survreg(fmla, data = data_fit, dist = dist)

  ## Fit the model for C with quantile_forest (now only supports 1d)
  fit_X <- data.frame(X = data_fit[,names(data_fit) %in% xnames])
  qc_mdl <- quantile_forest(fit_X, as.vector(data_fit$C))

  ## The fitted quantile for the calibration data
  xdf <- data.frame(data_calib[,names(data_calib) %in% xnames])
  colnames(xdf) = colnames(newdata)
  res <- predict(mdl,
                 newdata = xdf,
                 type = "quantile",
                 p = alpha)
  quant <-  res
  score <- quant-data_calib$censored_T


  ## Estimate the global censoring rate
  cens_rt <- 1 - sum(data_fit$event)/nrow(data_fit)

  cens_rt = 1-1/log(n)
  c_ref <- 1 : 6 / 2
  xmin <- 0; xmax <- 4
  exp_rate <- 0.4
  pr_all_list <- matrix(0, n + n_test, length(c_ref))
  mod <- "cox"
  num_grids <- 20
  len_x <- n_test

  ## Obtaining lower bounds
  simures = NULL
  simulen = NULL
  qt_res <- alpha_qt(mdl, newdata,
                     data_fit, data_calib,
                     xnames, alpha, len_x,
                     cens_rt = cens_rt,
                     num_grids = num_grids,
                     mdl0 = mdl0)
  res1 <- qt_res$lower_bnd_l

  output <- data.frame(qtl = res1)

  res0 = cfsurv(x,c_list=NULL,
                pr_list=NULL,
                pr_new_list=NULL,
                Xtrain,C,event,time,
                alpha=alpha,
                type="quantile",
                #seed = 24601,
                model = mod,
                dist= "weibull",
                I_fit = NULL,
                ftol=.1,tol=.1,
                n.tree=100)
  output$qc0 <- res0$res


  ########################################
  ### utility functinos
  ########################################
  ## A function to get result
  extract_res_univariate <- function(x){
    res <- mdl_coef[1] + mdl_coef[-1]*x
    return(res)
  }

  ## A utility function to extract quantiles from a coxph object
  extract_quant <- function(mdl, x, alpha){
    res <- summary(survfit(mdl, newdata = x))
    time_point <- res$time
    survcdf <- 1 - res$surv
    # if(sum(survcdf >= alpha)==0){
    #   quant = max(time_point)
    # }else{
    #   quant <- time_point[min(which(survcdf >= alpha))] 
    # }
    quant <- time_point[min(which(survcdf >= alpha))] 
    return(quant)
  }

  ########################################
  ## Quantile regression
  ########################################
  ## Cox
  cat("Compute the result of Cox...")
  fmla <- as.formula(paste("Surv(censored_T, event) ~ ",
                           paste(xnames, collapse= "+")))
  mdl <- coxph(fmla, data = data)
  res <- apply(data_test, 1, extract_quant, mdl = mdl,  alpha = alpha)
  output$cox.bnd <- res
  cat("Done.\n")


  cov_rt = function(lb_res,T_test){sum(T_test >= lb_res)/length(lb_res)}
  #covrt0 = apply(output, 2, cov_rt, T_test=T0)
  covrt = apply(output, 2, cov_rt, T_test=T_test)
  len = apply(output, 2, median)
  res = c(covrt, len, res0$c)#, covrt0)
  return(res)
}


library(foreach)
library(doParallel)
repN = 16
n = 500
alpha=0.1
detectCores()
cl<- makeCluster(16)
registerDoParallel(cl)
simures = foreach(i = 1:repN, .combine="rbind") %dopar% run.simu(i+1234,alpha,n)
stopCluster(cl)

names = c("Adaptive cutoff","Fixed cutoff","Cox")
colnames(simures) = c(rep(names,2),"c")

par(mfrow=c(1,2))
boxplot(simures[,1:3],col="white",border=2:4,
        ylab="Coverage rate",main="",ylim=c(0.8,1),
        cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
abline(h = 1-alpha, lty=2)
boxplot(simures[,4:6],col="white",border=2:4,
        ylab="Average LPB",main="",
        cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
par(mfrow=c(1,1))

sigma_x <- function(x) (x+5)/10
exp_rate = .1
gen_t <- function(x) rpareto(length(x), shape=4/x, location=0.1+beta*x)
gen_c <- function(x) rexp(rate = exp_rate, n = length(x))
# gen_c <- function(x) rpareto(length(x), shape=0.2, location=0.5)
#setwd("~/Documents/Github/adaptive_cfsurv-main/cfsuv/result")
#pdf("plot_1.pdf", width = 12, height = 4)
#par(mfrow = c(1,3))
n_test <- 100
n_train <- n
n_calib <- n
n = n_calib
beta <- 20 / sqrt(n)
xmin <- 0; xmax <- 4
########################################
## Generate training data
########################################
# set.seed(24601)
X <- runif(n_train, xmin, xmax)
T <- gen_t(X)
C <- gen_c(X)
print(mean(C<=T))
plot(T~X,col=2,pch=2,xlab="X",ylab="",#ylim=c(0,600),
     cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
lines(C~X,col=4,ty="p",pch=4)
legend("top", legend=c("log(survival time)","log(censoring time)"),
       pch=c(2,4), col = c(2,4), cex=1.3)
# names = c("Adaptive cutoff","Fixed cutoff","Cox")
# colnames(simures) = c(rep(names,2),"c")
# boxplot(simures[,1:3],col="white",border=2:4,
#         ylab="Coverage rate",main="",ylim=c(0.7,1),
#         cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
# abline(h = 1-alpha, lty=2)
# boxplot(simures[,4:6],col="white",border=2:4,
#         ylab="Average LPB",main="",ylim = c(0,12),
#         cex.lab=1.3, cex.axis=1.3, cex.main=1.3, cex.sub=1.3)
# par(mfrow = c(1,3))
#dev.off()

