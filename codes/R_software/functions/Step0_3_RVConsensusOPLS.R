#' RVConsensusOPLS
#' Consensus OPLS-DA with RV coefficients weighting and 
#' DQ2 computation for discriminant analysis
#' 
#' @param data: the collection list containing each block of data.
#' @param Y: The response matrix (un-centered/scaled).
#' @param A: The number of Y-predictive components (integer). 
#' @param maxOrtholvs: The maximal number of Y-orthogonal components (integer).
#' @param nrcv: Number of cross-validation rounds (integer).
#' @param cvType: Type of cross-validation used. Either `nfold` for n-fold
#' cross-validation, `mccv` for Monte Carlo CV or `mccvb` for Monte Carlo 
#' class-balanced CV.
#' @param modelType: type of OPLS regression model. Can be defined as "reg" for 
#' regression or "da" for discriminant analysis. Default value "da".
#' @param verbose: logical which indicates whether the user wants to see the 
#' progress bar printlayed in the ConsensusOLPSCV function.
#'
#' @return 
#' `execution_time`: function execution time.
#' `model`: a list with all model parameters.
#'
#' @examples
#' TO DO

RVConsensusOPLS <- function(data = collection,
                            Y,
                            A = 1, 
                            maxOrtholvs = 10, 
                            nrcv = 100,
                            cvType = "nfold",
                            modelType = "da", 
                            verbose = FALSE){
  # Variable format control
  if(!is.list(data)){stop("data is not a list.")}
  if(!is.matrix(Y)){stop("Y is not a matrix.")}
  if(!is.numeric(A)){stop("A is not numeric.")}
  if(!is.numeric(maxOrtholvs)){stop("maxOrtholvs is not numeric.")}
  if(!is.numeric(nrcv)){stop("nrcv is not numeric.")}
  if(!is.character(cvType)){
    stop("cvType is not a character.")
  } else{
    if(!(cvType %in% c("nfold", "mccv", "mccvb"))){
      stop("cvType must be `nfold`, `mccv` or `mccvb`.")
    }
  }
  if(!is.character(modelType)){
    stop("modelType is not a character.")
  } else{
    if(!(modelType %in% c("reg", "da"))){
      stop("modelType must be `reg` or `da`.")
    }
  }
  if(!is.logical(verbose)){stop("verbose must be logical `TRUE` or `FALSE`.")}
  
  # Function loading control
  if (!exists("koplsScale", mode = "function")) {
    warning("Remember to load the source code for the `koplsScale` function.")
  }
  if (!exists("koplsKernel", mode = "function")) {
    warning("Remember to load the source code for the `koplsKernel` function.")
  }
  if (!exists("RVmodified", mode = "function")) {
    warning("Remember to load the source code for the `RVmodified` function.")
  }
  if (!exists("ConsensusOPLSCV", mode = "function")) {
    warning("Remember to load the source code for the `ConsensusOPLSCV` function.")
  }
  
  # Evaluate time ellapse
  tStart <- Sys.time()
  
  # Check collection dimension
  ntable <- base::length(collection)
  nrow <- nrow(collection[[1]])
  
  # Initialize parameters
  W_mat <- base::matrix(data = 0, nrow = nrow, ncol = nrow)
  preProcK <- "mc"
  preProcY <- "mc"
  
  #Fraction of data to be used for cross-validation
  cvFrac <- 0.75
  
  if (modelType == "reg") {
    koplsScale <- koplsScale(X = Y, centerType = preProcY, scaleType = "no")
    Yc <- koplsScale$matrix
  } else {
    Yc <- Y
  }
  
  # For each data block
  xnorm <- list() ; AMat <- list() ; RV <- list()
  for (ta in 1:ntable) {
    # Produce the kernel of the data block
    temp <- koplsKernel(X1 = collection[[ta]], X2 = NULL, Ktype = 'p', params = 1)
    # Frobenius norm of the kernel
    xnorm[[ta]] <- base::norm(x = temp, type = "F")
    # Normalize the Kernel
    AMat[[ta]] <- temp/xnorm[[ta]]
    # RV coefficient for AMat
    RV[[ta]] <- (RVmodified(X = AMat[[ta]], Y = Yc) + 1) / 2
    # calculates the weighted sum of blocks kernel by the RV coeff
    W_mat <- W_mat + RV[[ta]] * AMat[[ta]]
  }
  
  # Performs a Kernel-OPLS cross-validation for W_mat
  modelCV <- ConsensusOPLSCV(K = W_mat, Y = Y, A = A, oax = maxOrtholvs, 
                             nbrcv = nrcv, cvType = cvType, preProcK = preProcK, 
                             preProcY = preProcY, cvFrac = cvFrac, 
                             modelType = modelType, verbose = verbose)
  
  Ylarg <- ncol(Y)  
  
  # Search for the optimal model based on DQ2
  if (modelType == 'da') {
    dqq <- base::matrix(data = 0, nrow = maxOrtholvs+1, ncol = Ylarg)
    PRESSD <- base::matrix(data = 0, nrow = maxOrtholvs+1, ncol = Ylarg)
    
    for (i in 0:maxOrtholvs) {
      for (j in 1:Ylarg) {
        # For each Y column, perform the DQ2
        result <- DQ2(Ypred = base::matrix(data = modelCV$cv$AllYhat[, Ylarg*i+j],
                                           ncol = 1), 
                      Y = base::matrix(data = Y[, j], 
                                       ncol = 1))
        dqq[i+1, j] <- result$dqq  
        PRESSD[i+1, j] <- result$PRESSD
      }
    }
    
    dq2 <- base::apply(X = dqq, MARGIN = 1, FUN = function(X) mean(X))
    index <- A  
    
    # Finds the optimal number of orthogonal components as a function of DQ2
    while (index < (maxOrtholvs+A) && (dq2[index+1] - dq2[index]) > 0.01) {
      index <- index + 1
    }
    
    # Add DQ2 in the model objects
    modelCV$cv$DQ2Yhat <- dq2 
    # Add optimal number of orthogonal components in the model objects
    modelCV$cv$OrthoLVsOptimalNum <- index - A
    
  } else { # if modelType == "reg"
    index <- A 
    
    # Finds the optimal number of orthogonal components as a function of Q2Yhat
    while (index < (maxOrtholvs+A) && 
           (modelCV$cv$Q2Yhat[index+1] - modelCV$cv$Q2Yhat[index]) > 0.01) {
      index <- index + 1
    }
    # Add optimal number of orthogonal components in the model objects
    modelCV$cv$OrthoLVsOptimalNum <- index - A
  }
  
  # Simplifies the name to be used afterwards
  if (modelCV$cv$OrthoLVsOptimalNum == 0) {
    OrthoLVsNum <- 1
  } else {
    OrthoLVsNum <- modelCV$cv$OrthoLVsOptimalNum
  }
  
  # Recompute the optimal model using OrthoLVsNum parameters
  modelCV$Model <- koplsModel(K = W_mat, Y = Y, A = A, nox = OrthoLVsNum, 
                              preProcK = preProcK, preProcY = preProcY)
  
  # Adjust Yhat to the selected model size
  modelCV$cv$Yhat <- modelCV$cv$AllYhat[, ((Ylarg*A)+(OrthoLVsNum*A)) : 
                                          ((Ylarg*A)+(OrthoLVsNum*A)+Ylarg-1)]
  
  # Compute the blocks contributions for the selected model
  lambda <- base::matrix(data = 0, nrow = ntable, ncol = A+OrthoLVsNum)
  for (j in 1:ntable) {
    for (k in 1:A) {
      T <- modelCV$Model$T[, k]
      lambda[j, k] <- t(T) %*% AMat[[j]] %*% T
    }
    for (l in 1:OrthoLVsNum) {
      To <- modelCV$Model$To[, l]
      lambda[j, l+A] <- t(To) %*% AMat[[j]] %*%To
    }
  }
  # Stores raw lambda coefficient values in the model object
  modelCV$Model$lambda_raw <- lambda 
  
  # Normalize the lambda coefficients
  for (nb in 1:ncol(lambda)) {
    lambda[, nb] <- lambda[, nb] / base::sum(lambda[, nb])
  }
  # Stores normalized lambda values in the model object
  modelCV$Model$lambda <- lambda 
  
  # Compute the loadings for the selected model size
  loadings <-  base::matrix(data = list(), nrow = ntable, ncol = (A+OrthoLVsNum))
  for (ta in 1:ntable) {
    for (m in 1:A) {
      T <- base::matrix(modelCV$Model$T[, m], ncol = 1)
      loadings[[ta, m]] <- t(collection[[ta]]) %*% T %*% base::solve(t(T) %*% T)
    }
    for (n in 1:OrthoLVsNum) {
      To <- modelCV$Model$To[, n]
      loadings[[ta, n+m]] <- list(t(collection[[ta]]) %*% To %*% base::solve(t(To) %*% To))
    }
  }
  
  # Add RV coefficients in the model objects
  modelCV$RV <- RV  
  # Add the loadings in the model objects
  modelCV$Model$loadings <- loadings 
  
  tStop <- Sys.time()
  return(list("execution_time" = as.numeric(tStop - tStart), 
              "model" = modelCV))
}