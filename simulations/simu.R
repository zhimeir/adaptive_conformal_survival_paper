simu <- function(seed, setting, n, p, 
                 n_train, n_calib, n_test, 
                 beta, xmin, xmax, 
                 exp_rate, alpha, mod_list){
  
  

  ## Initialization
  set.seed(seed)
  xnames <- paste0("X",1:p) 

  ## Generate data according to the setting
  data_obj <- model_generating_fun(n_train, n_calib, n_test,
                                   setting, beta, xnames, 
                                   xmin, xmax, exp_rate) 
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
  fit = data_fit
  fit$C = -fit$C
  alpha_list = alpha

  ## Obtaining lower bounds
  simures = NULL
  simulen = NULL
  for(mod in mod_list){
    cat("Compute the result of qt...\n")
    cat("Compute the result of qct...\n")
    start_time = proc.time()[3]
    mdl0 = GauPro(X = as.matrix(fit[,names(fit) %in% xnames]),
                  Z = fit$C, D = p,
                  type = "Gauss")
    lb_res = cfsurv_c(x=x,
                      Xtrain = Xtrain,
                      C = C,
                      event = event,
                      time=time,
                      alpha=alpha,
                      mdl0=mdl0)
    end_time <- proc.time()[3]
    time_q <- end_time - start_time
    
    cat("Compute the result of qc0...\n")
    start_time <- proc.time()[3]
    res0 = cfsurv(x,c_list=NULL,
                  pr_list=NULL,
                  pr_new_list=NULL,
                  Xtrain,C,event,time,
                  alpha=alpha,
                  type="quantile",
                  model = mod,
                  dist= "weibull",
                  I_fit = NULL,
                  ftol=.1,tol=.1,
                  n.tree=100)
    end_time <- proc.time()[3]
    time_qc0 <- end_time - start_time
    
    cov_rt = function(lb_res,T_test){sum(T_test >= lb_res)/length(lb_res)}
    
    res1 = lb_res$lower_bnd_qtl
    res2 = lb_res$lower_bnd_qctl
    time_qt = lb_res$time_qt
    time_qct = lb_res$time_qct
    time_q_base = time_q - (time_qt+time_qct)
    time_qt = time_qt + time_q_base
    time_qct = time_qct + time_q_base
    
    time <- c(time_qt)
    time <- c(time, time_qct)
    time <- c(time, time_qc0)
    output <- data.frame(qtl = res1)
    output$qctl <- res2
    output$qc0 <- res0$res
    
    # parametric methods
    ########################################
    ## vanilla CQR
    ########################################
    cat(" - Computing the result of vanilla-CQR...\n")
    start_time <- proc.time()[3]
    res <- lapply(alpha_list,
                  cqr,
                  x=data_test[,names(data_test) %in% xnames],
                  Xtrain = data[,names(data) %in% xnames],
                  Ytrain = data$censored_T,
                  I_fit = 1:n_train,
                  seed = seed +7)
    res <- do.call(rbind,lapply(res,as.data.frame))
    output$cqr.bnd <- res[,1]
    end_time <- proc.time()[3]
    time <- c(time, end_time - start_time)
    cat("done.\n")
    
    
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
      if(sum(survcdf >= alpha)==0){
        quant = max(time_point)
      }else{
        quant <- time_point[min(which(survcdf >= alpha))]
      }
      # quant <- time_point[min(which(survcdf >= alpha))]
      return(quant)
    }
    
    ########################################
    ## Quantile regression
    ########################################
    ## Cox
    cat("Compute the result of Cox...\n")
    start_time <- proc.time()[3]
    fmla <- as.formula(paste("Surv(censored_T, event) ~ ",
                             paste(xnames, collapse= "+")))
    mdl <- coxph(fmla, data = data)
    res <- c()
    for(alpha in alpha_list){
      res <- c(res, apply(data_test, 1, extract_quant, mdl = mdl,  alpha = alpha))
    }
    output$cox.bnd <- res
    end_time <- proc.time()[3]
    time <- c(time, end_time - start_time)
    cat("Done.\n")
    
    
    ## random forest
    cat("Compute the result of random forest...")
    start_time <- proc.time()[3]
    ntree <- 1000
    nodesize <- 80
    res <- c()
    fmla <- as.formula(paste("censored_T ~ ",
                             paste(xnames, collapse= "+")))
    for(alpha in alpha_list){
      if(p==1){
        mdl <- crf.km(as.formula("censorerd_T~X1"),
                      ntree = ntree, 
                      nodesize = nodesize,
                      data_train = data[,colnames(data) %in% c("X1","censored_T","event")], 
                      data_test = data.frame(X1 = data_test$X1), 
                      yname = 'censored_T', 
                      iname = 'event',
                      tau = alpha,
                      method = "grf")
      }else{
        mdl <- crf.km(fmla,
                      ntree = ntree, 
                      nodesize = nodesize,
                      data_train = data[,colnames(data) %in% c(xnames,"censored_T","event")], 
                      data_test = data_test[,colnames(data_test) %in% xnames], 
                      yname = 'censored_T', 
                      iname = 'event',
                      tau = alpha,
                      method = "grf")
      }
      res <- c(res,mdl$predicted)
    }
    output$rf.bnd  <- res
    end_time <- proc.time()[3]
    time <- c(time, end_time - start_time)
    cat("Done.\n")
    
    simures = apply(output, 2, cov_rt, T_test=T_test)
    simulen = apply(output, 2, mean)
  }
  simu_out = c(simures,simulen,time)
  return(simu_out=simu_out)
}


