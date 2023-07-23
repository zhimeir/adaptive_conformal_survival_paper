# function to construct conformal confidence interval
cfsurv_c <- function(x,
                     Xtrain, C, event, time,
                     alpha=0.05,
                     I_fit = NULL,mdl0
){

  ## Process the input
  ## Check the length of x and c: only two cases are supported. length(r)=1, or length(r)=length(x)
  X <- Xtrain
  if(is.null(dim(x)[1])){
    len_x <- length(x)
    p <- 1
  }else{
    len_x <- dim(x)[1]
    p <- dim(x)[2]
  }
  
  
  if(is.null(dim(X)[1])){
    n <- length(X)
    pX <- 1
  }else{
    n <- dim(X)[1]
    pX <- dim(X)[2]
  }

  ## Check the dimensions of the data 
  xnames <- paste0('X', 1:p)
  data <- as.data.frame(cbind(C,event,time,X))
  colnames(data) <- c("C","event","censored_T",xnames)
  
  ## Split the data into the training set and the calibration set
  n <- dim(data)[1]
  n_train <- n/2
  n_calib <- n - n_train
  if(is.null(I_fit)){
    I_fit <- sample(1:n, n_train, replace = FALSE)
  }
  data_fit <- data[I_fit,]
  data_calib <- data[-I_fit,]
  
  ## Run the main function with Cox base algorithm
  res <- cox_based(x, alpha,
                   data_fit,
                   data_calib,
                   dist, mdl0)
  

  return(res)
  
  
}

