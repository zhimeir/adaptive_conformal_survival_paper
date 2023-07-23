############################################
## estimate miscoverage rate 
## using estimated quantile of T
############################################
alpha_qt <- function(mdl, newdata,
                     data_fit, data_calib,
                     xnames, alpha, len_x,
                     mdl0, cens_rt){
  
  v_list = v_pts_qt(mdl,
                      data_fit, data_calib,
                      xnames, alpha, cens_rt)
  v_list = sort(unique(as.numeric(v_list)))
  
  ## Obtain the final confidence interval
  lower_bnd_l <- rep(0,len_x)
  lower_bnd_g <- rep(0,len_x)
  
  alpha_v_list <- sapply(v_list, est_alpha_qt, mdl = mdl, 
                    data_calib = data_calib, xnames = xnames, alpha = alpha, 
                    cens_rt = cens_rt, mdl0 = mdl0, newdata = newdata)

  # monotonize alpha
  alpha_v <- monot(alpha_v_list)
  # return 0 if alpha is above target level
  if(sum(alpha_v<=alpha)==0){
    v_hat_l = NULL
    v_hat_g = NULL
  }else{
    v_hat_l <- min(v_list[alpha_v <= alpha])
  }
  
  for(i in 1:len_x){
    nxi <- data.frame(newdata[i,])
    colnames(nxi) <- xnames
    
    lower_bnd_l[i] <- lv_qt(mdl, nxi, v_hat_l, alpha, cens_rt)  
  }
  
  return(list(lower_bnd_l = lower_bnd_l, 
              lower_bnd_g = lower_bnd_g))
  
}

############################################
## compute estimated miscoverage rate
############################################
est_alpha_qt <- function(mdl, data_calib, xnames, v, alpha, cens_rt, mdl0, newdata){
  ## Check the dimensionality of the input
  if(is.null(dim(newdata)[1])){
    len_x <- length(newdata)
    p <- 1
  }else{
    len_x <- dim(newdata)[1]
    p <- dim(newdata)[2]
  }
  calib_x <- data.frame(X = data_calib[,names(data_calib) %in% xnames])
  names(calib_x) = xnames
  
  lv_calib = lv_qt(mdl, calib_x, v, alpha, cens_rt)
  cens = cens_prob(mdl0,data_calib,NULL,
                   method = "gpr",
                   xnames = xnames,
                   c = lv_calib)
  pr_calib = cens$pr_calib
  weight_calib <- 1/pr_calib
  
  w_new = 0
  
  ind1 = (data_calib$censored_T < lv_calib) & (data_calib$C >= lv_calib)
  ind2 = (data_calib$C >= lv_calib)
  
  sum_num <- sum(weight_calib[ind1]) + w_new
  sum_den <- sum(weight_calib[ind2]) + w_new
  if((length(ind2)==0)||is.na(sum_num/sum_den)){
    alphav = 1
  }else{
    alphav <- sum_num / sum_den
  }
  
  return(alphav)
} 

############################################
## compute lower prediction bound
############################################
lv_qt <- function(mdl, calib_x, v, alpha, cens_rt){
  if(length(v)==0){
    return(lv2_calib = 0)
  }
  lv2_calib <- predict(mdl, newdata = calib_x, type = "quantile", p = 1-v)
  return(lv2_calib)
}


################################################################
## The function returns  predictive intervals resulting
## from the L_v defined based on integrtaed quantiles
################################################################
alpha_qct <- function(mdl, qc_mdl, newdata,
                     data_fit, data_calib,
                     xnames, alpha, len_x,
                     mdl0, cens_rt){

  v_list = v_pts_qct(mdl, qc_mdl,
                      data_fit, data_calib,
                      xnames, alpha, cens_rt)
  v_list = sort(unique(as.numeric(v_list)))

  ## Obtain the final confidence interval
  lower_bnd_l <- rep(0,len_x)
  lower_bnd_g <- rep(0,len_x)
  
  alpha_v_list <- sapply(v_list, est_alpha_qct, mdl = mdl, qc_mdl = qc_mdl,
                    data_calib = data_calib, xnames = xnames, alpha = alpha, 
                    cens_rt = cens_rt, mdl0 = mdl0, newdata = newdata)
  
  # monotonize alpha
  alpha_v <- monot(alpha_v_list)
  # return 0 if alpha is above target level
  if(sum(alpha_v<=alpha)==0){
    v_hat_l = NULL
    v_hat_g = NULL
  }else{
    v_hat_l <- min(v_list[alpha_v <= alpha])
  }
  
  for(i in 1:len_x){
    nxi <- data.frame(newdata[i,])
    colnames(nxi) <- xnames

    lower_bnd_l[i] <- lv_qct(mdl, qc_mdl, nxi, v_hat_l, alpha, cens_rt)
  }

  return(list(lower_bnd_l = lower_bnd_l, 
              lower_bnd_g = lower_bnd_g))

}


################################################################
## Compute estimated alpha
################################################################
est_alpha_qc <- function(mdl, data_calib, xnames, v){

  calib_x <- data.frame(X = data_calib[,names(data_calib) %in% xnames])
  names(calib_x) = xnames

  lv_calib <- predict(mdl, calib_x, 1 - v)
  lv_calib <- (lv_calib$predictions)[,1]

  sum_num <- sum((data_calib$censored_T < lv_calib) * (data_calib$C >= lv_calib)) + 1
  sum_den <- sum((data_calib$C >= lv_calib)) + 1
  
  alphav <- sum_num / sum_den

  return(alphav)
} 


################################################################
## Compute estimated alpha using integrated quantiles
################################################################
est_alpha_qct <- function(mdl, qc_mdl, data_calib, xnames, v, alpha, cens_rt, mdl0, newdata){
  ## Check the dimensionality of the input
  if(is.null(dim(newdata)[1])){
    len_x <- length(newdata)
    p <- 1
  }else{
    len_x <- dim(newdata)[1]
    p <- dim(newdata)[2]
  }
  calib_x <- data.frame(X = data_calib[,names(data_calib) %in% xnames])
  names(calib_x) = xnames

  lv_calib = lv_qct(mdl, qc_mdl, calib_x, v, alpha, cens_rt)
  
  cens = cens_prob(mdl0,data_calib,NULL,
                   method = "gpr",
                   xnames = xnames,
                   c = lv_calib)
  pr_calib = cens$pr_calib
  weight_calib <- 1/pr_calib
  
  w_new = 0
  
  ind1 = (data_calib$censored_T < lv_calib) & (data_calib$C >= lv_calib)
  ind2 = (data_calib$C >= lv_calib)
  
  sum_num <- sum(weight_calib[ind1]) + w_new
  sum_den <- sum(weight_calib[ind2]) + w_new
  if((length(ind2)==0)||is.na(sum_num/sum_den)){
    alphav = 1
  }else{
    alphav <- sum_num / sum_den
  }

  return(alphav)
} 

lv_qct <- function(mdl, qc_mdl, calib_x, v, alpha, cens_rt){
  if(length(v)==0){
    return(lv_calib = 0)
  }
  lv1_calib <- predict(qc_mdl, calib_x, cens_rt)
  lv2_calib <- predict(mdl, newdata = calib_x, type = "quantile", p = 1-v)
  lv1_calib <- (lv1_calib$predictions)[,1]
  lv_calib = pmin(lv1_calib, lv2_calib)
  return(lv_calib)
}


############################################
## determine the finite candidate set for v
############################################
inverse <- function(f, lower, upper){
  function(y){
    uniroot(function(x){f(x) - y}, lower = lower, upper = upper, tol=1e-5)[1]
  }
}

v_pts_qt = function(mdl,
                      data_fit, data_calib,
                      xnames, alpha, cens_rt){
  pts = matrix(0,nrow(data_calib),2)
  for(i_calib in 1:nrow(data_calib)){
    calib_x = data.frame(X = data_calib[i_calib,names(data_calib) %in% xnames])
    names(calib_x) = xnames
    lv = function(v) lv_qt(mdl, calib_x, v, alpha, cens_rt)
    x0 = 0.02
    x1 = 0.95
    upb = lv(x0)
    lwb = lv(x1)
    lv_inv = inverse(function(v) lv_qt(mdl, calib_x, v, alpha, cens_rt), x0, x1)
    if(data_calib$event[i_calib]==0){
      if(data_calib$C[i_calib] >= upb){
        pts[i_calib,1] = pts[i_calib,2] = x0
      }
      if(data_calib$C[i_calib] <= lwb){
        pts[i_calib,1] = pts[i_calib,2] = x1
      }
      if((data_calib$C[i_calib] > lwb) && (data_calib$C[i_calib] < upb)){
        pts[i_calib,1] = pts[i_calib,2] = lv_inv(data_calib$C[i_calib])$root
      }
    }
    if(data_calib$event[i_calib]==1){
      if(data_calib$C[i_calib] >= upb){
        pts[i_calib,1] = x0
      }
      if(data_calib$C[i_calib] <= lwb){
        pts[i_calib,1] = x1
      }
      if((data_calib$C[i_calib] > lwb) && (data_calib$C[i_calib] < upb)){
        pts[i_calib,1] = lv_inv(data_calib$C[i_calib])$root
      }
      if(data_calib$censored_T[i_calib] >= upb){
        pts[i_calib,2] = x0
      }
      if(data_calib$censored_T[i_calib] <= lwb){
        pts[i_calib,2] = x1
      }
      if((data_calib$censored_T[i_calib] > lwb) && (data_calib$censored_T[i_calib] < upb)){
        pts[i_calib,2] = lv_inv(data_calib$censored_T[i_calib])$root
      }
    }
  }
  return(pts)
}

v_pts_qc = function(qc_mdl,
                    data_fit, data_calib,
                    xnames, alpha, cens_rt){
  pts = matrix(0,nrow(data_calib),2)
  for(i_calib in 1:nrow(data_calib)){
    calib_x = data.frame(X = data_calib[i_calib,names(data_calib) %in% xnames])
    names(calib_x) = xnames
    lv = function(v) lv_qc(qc_mdl, calib_x, v)
    x0 = 0.02
    x1 = 0.95
    upb = lv(x0)
    lwb = lv(x1)
    lv_inv = inverse(function(v) lv_qc(qc_mdl, calib_x, v), x0, x1)
    if(data_calib$event[i_calib]==0){
      if(data_calib$C[i_calib] >= upb){
        pts[i_calib,1] = pts[i_calib,2] = x0
      }
      if(data_calib$C[i_calib] <= lwb){
        pts[i_calib,1] = pts[i_calib,2] = x1
      }
      if((data_calib$C[i_calib] > lwb) && (data_calib$C[i_calib] < upb)){
        pts[i_calib,1] = pts[i_calib,2] = lv_inv(data_calib$C[i_calib])$root
      }
    }
    if(data_calib$event[i_calib]==1){
      if(data_calib$C[i_calib] >= upb){
        pts[i_calib,1] = x0
      }
      if(data_calib$C[i_calib] <= lwb){
        pts[i_calib,1] = x1
      }
      if((data_calib$C[i_calib] > lwb) && (data_calib$C[i_calib] < upb)){
        pts[i_calib,1] = lv_inv(data_calib$C[i_calib])$root
      }
      if(data_calib$censored_T[i_calib] >= upb){
        pts[i_calib,2] = x0
      }
      if(data_calib$censored_T[i_calib] <= lwb){
        pts[i_calib,2] = x1
      }
      if((data_calib$censored_T[i_calib] > lwb) && (data_calib$censored_T[i_calib] < upb)){
        pts[i_calib,2] = lv_inv(data_calib$censored_T[i_calib])$root
      }
    }
  }
  return(pts)
}

v_pts_qct = function(mdl, qc_mdl,
                    data_fit, data_calib,
                    xnames, alpha, cens_rt){
  pts = matrix(0,nrow(data_calib),2)
  for(i_calib in 1:nrow(data_calib)){
    calib_x = data.frame(X = data_calib[i_calib,names(data_calib) %in% xnames])
    names(calib_x) = xnames
    lv = function(v) lv_qct(mdl, qc_mdl, calib_x, v, alpha, cens_rt)
    x0 = 0.02
    x1 = 0.95
    upb = lv(x0)
    lwb = lv(x1)
    lv_inv = inverse(function(v) lv_qct(mdl, qc_mdl, calib_x, v, alpha, cens_rt), x0, x1)
    if(data_calib$event[i_calib]==0){
      if(data_calib$C[i_calib] >= upb){
        pts[i_calib,1] = pts[i_calib,2] = x0
      }
      if(data_calib$C[i_calib] <= lwb){
        pts[i_calib,1] = pts[i_calib,2] = x1
      }
      if((data_calib$C[i_calib] > lwb) && (data_calib$C[i_calib] < upb)){
        pts[i_calib,1] = pts[i_calib,2] = lv_inv(data_calib$C[i_calib])$root
      }
    }
    if(data_calib$event[i_calib]==1){
      if(data_calib$C[i_calib] >= upb){
        pts[i_calib,1] = x0
      }
      if(data_calib$C[i_calib] <= lwb){
        pts[i_calib,1] = x1
      }
      if((data_calib$C[i_calib] > lwb) && (data_calib$C[i_calib] < upb)){
        pts[i_calib,1] = lv_inv(data_calib$C[i_calib])$root
      }
      if(data_calib$censored_T[i_calib] >= upb){
        pts[i_calib,2] = x0
      }
      if(data_calib$censored_T[i_calib] <= lwb){
        pts[i_calib,2] = x1
      }
      if((data_calib$censored_T[i_calib] > lwb) && (data_calib$censored_T[i_calib] < upb)){
        pts[i_calib,2] = lv_inv(data_calib$censored_T[i_calib])$root
      }
    }
  }
  return(pts)
}
