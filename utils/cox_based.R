############################################
## Lower prediction bound based on Cox model
## with adaptive cutoffs
############################################

cox_based <- function(x,alpha,
                      data_fit,
                      data_calib,
                      dist, mdl0){
  
  ## Check the dimensionality of the input
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }
  n <- nrow(data_calib)
  ## Fit the survival model
  xnames <- paste0("X",1:p)
  newdata <- data.frame(x)
  colnames(newdata) <- xnames
  fmla <- as.formula(paste("Surv(censored_T, event) ~ ", paste(xnames, collapse= "+")))
  mdl <- survreg(fmla, data = data_fit, dist= "weibull")

  ## The fitted quantile for the calibration data
  xdf <- data.frame(data_calib[,names(data_calib) %in% xnames])
  colnames(xdf) = colnames(newdata)
  res <- predict(mdl,
                 newdata = xdf,
                 type = "quantile",
                 p = alpha)
  quant <-  res  
  score <- quant-data_calib$censored_T
  
  # cutoff for quantile of C
  # to bound weights
  cens_rt = 1-1/log(n)
  
  start_time = proc.time()[3]
  qt_res <- alpha_qt(mdl, newdata,
                     data_fit, data_calib,
                     xnames, alpha, len_x,
                     cens_rt = cens_rt,
                     mdl0 = mdl0)
  lower_bnd_qtg <- qt_res$lower_bnd_g
  lower_bnd_qtl <- qt_res$lower_bnd_l
  end_time <- proc.time()[3]
  time_qt <- end_time - start_time
  
  start_time = proc.time()[3]
  ## Fit the model for C with quantile_forest (now only supports 1d)
  fit_X <- data.frame(X = data_fit[,names(data_fit) %in% xnames])
  qc_mdl <- quantile_forest(fit_X, as.vector(data_fit$C))
  
  qct_res <- alpha_qct(mdl, qc_mdl, newdata,
                     data_fit, data_calib,
                     xnames, alpha, len_x,
                     cens_rt = cens_rt,
                     mdl0 = mdl0)
  lower_bnd_qctg <- qct_res$lower_bnd_g
  lower_bnd_qctl <- qct_res$lower_bnd_l
  end_time <- proc.time()[3]
  time_qct <- end_time - start_time
  

  return(list(lower_bnd_qtg =  lower_bnd_qtg,
              lower_bnd_qtl = lower_bnd_qtl,
              lower_bnd_qctg = lower_bnd_qctg,
              lower_bnd_qctl = lower_bnd_qctl,
              time_qt = time_qt,
              time_qct = time_qct
              ))
}



