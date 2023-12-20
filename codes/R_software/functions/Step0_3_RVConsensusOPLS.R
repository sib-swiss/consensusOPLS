#' @title RVConsensusOPLS
#' @description
#' Consensus OPLS-DA with RV coefficients weighting and 
#' DQ2 computation for discriminant analysis.
#' 
#' @param data the collection list containing each block of data.
#' @param Y The response matrix (un-centered/scaled).
#' @param A The number of Y-predictive components (integer). 
#' @param maxOrtholvs The maximal number of Y-orthogonal components (integer).
#' @param nrcv Number of cross-validation rounds (integer).
#' @param cvType Type of cross-validation used. Either \code{nfold} for n-fold
#' cross-validation, \code{mccv} for Monte Carlo CV or \code{mccvb} for Monte 
#' Carlo class-balanced CV.
#' @param modelType type of OPLS regression model. Can be defined as \code{reg} 
#' for regression or \code{da} for discriminant analysis. Default value is
#' \code{da}.
#' @param verbose logical which indicates whether the user wants to see the 
#' progress bar printed in the \code{ConsensusOLPSCV} function.
#'
#' @return A consensus OPLS model
#'
#' @examples
#' data(demo_3_Omics)
#' RVConsensusOPLS(data = demo_3_Omics[c("MetaboData", "MicroData","ProteoData")], 
#'                 Y = demo_3_Omics$Y[,1,drop=F])
#' 
#' @keywords internal   

RVConsensusOPLS <- function(data,
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

  # Check collection dimension
  ntable <- base::length(data)
  nrow <- nrow(data[[1]])
  
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
  RA <- mclapply(X = 1:ntable, mc.cores = detectCores(), FUN = function(i) {
    # Produce the kernel of the data block
    tmp <- koplsKernel(X1 = data[[i]], X2 = NULL, Ktype = 'p', params = c(order=1))
    # Frobenius norm of the kernel
    xnorm <- norm(x=tmp, type='F')
    # Normalize the Kernel
    AMat <- tmp/xnorm
    # RV coefficient for AMat
    RV <- (RVmodified(X = AMat, Y = Yc) + 1) / 2
    
    return (list(RV=RV, AMat=AMat))
  })
  # RV coefficient 
  RV <- mclapply(RA, mc.cores = detectCores(), FUN = function(x) x$RV)
  # Normalized kernel
  AMat <- mclapply(RA, mc.cores = detectCores(), FUN = function(x) x$AMat)
  # calculates the weighted sum of blocks kernel by the RV coeff
  W_mat <- Reduce("+", mclapply(1:ntable, mc.cores = detectCores(), FUN = function(i) {
    RA[[i]]$RV * RA[[i]]$AMat
  }))
  
  # Performs a Kernel-OPLS cross-validation for W_mat
  modelCV <- ConsensusOPLSCV(K = W_mat, Y = Y, A = A, oax = maxOrtholvs, 
                             nbrcv = nrcv, cvType = cvType, preProcK = preProcK, 
                             preProcY = preProcY, cvFrac = cvFrac, 
                             modelType = modelType, verbose = verbose)
  
  Ylarg <- ncol(Y)  
  
  # Search for the optimal model based on DQ2
  if (modelType == 'da') {
    mc <- (maxOrtholvs+1)*Ylarg
    mcj <- min(sqrt(mc), Ylarg)
    mci <- floor(detectCores()/mcj)
    results <- mclapply(X = 0:maxOrtholvs, mc.cores = mci, FUN = function(i) {
      mclapply(X = 1:Ylarg, mc.cores = mcj, FUN = function(j) {
        # For each Y column, perform the DQ2
        result <- DQ2(Ypred = matrix(data = modelCV$cv$AllYhat[, Ylarg*i+j],
                                     ncol = 1), 
                      Y = Y[, j, drop=F])
        return (result)
        
      })
    })
    dqq <- do.call(rbind, mclapply(X = 0:maxOrtholvs, mc.cores = mci, FUN = function(i) {
      do.call(cbind, mclapply(X = 1:Ylarg, mc.cores = mcj, FUN = function(j) {
        return (results[[i]][[j]]$dqq)
      }))
    }))
    PRESSD <- do.call(rbind, mclapply(X = 0:maxOrtholvs, mc.cores = mci, FUN = function(i) {
      do.call(cbind, mclapply(X = 1:Ylarg, mc.cores = mcj, FUN = function(j) {
        return (results[[i]][[j]]$PRESSD)
      }))
    }))
    
    dq2 <- rowMeans(dqq)
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
  modelCV$cv$Yhat <- modelCV$cv$AllYhat[, ((Ylarg*A)+(OrthoLVsNum*A)) + 0:(Ylarg-1), drop=F]
  
  # Compute the blocks contributions for the selected model
  lambda <- cbind(do.call(rbind,
                          mclapply(X = 1:ntable, mc.cores = detectCores(), FUN = function(j) {
                            diag(crossprod(modelCV$Model$T[, 1:A], 
                                           crossprod(t(AMat[[j]]), modelCV$Model$T[, 1:A])))
                          })),
                  do.call(rbind,
                          mclapply(X = 1:ntable, mc.cores = detectCores(), FUN = function(j) {
                            diag(crossprod(modelCV$Model$To[, 1:OrthoLVsNum], 
                                           crossprod(t(AMat[[j]]), modelCV$Model$To[, 1:OrthoLVsNum])))
                          })))
  
  # Stores raw lambda coefficient values in the model object
  modelCV$Model$lambda_raw <- lambda 
  
  # Normalize the lambda coefficients
  lambda <- sweep(x = lambda, MARGIN = 2, STATS = colSums(lambda), FUN = '/')
  
  # Stores normalized lambda values in the model object
  modelCV$Model$lambda <- lambda 
  
  # Compute the loadings for the selected model size
  loadings.list <- list( 
    mclapply(X = 1:ntable, mc.cores = detectCores(), FUN = function(i) {
      tcrossprod(crossprod(data[[i]], 
                           modelCV$Model$T[, 1:A, drop=F]), 
                 diag(diag(crossprod(modelCV$Model$T[, 1:A, drop=F]))^-1, ncol=A))
    }),
    mclapply(X = 1:ntable, mc.cores = detectCores(), FUN = function(i) {
      tcrossprod(crossprod(data[[i]],
                           modelCV$Model$To[, 1:OrthoLVsNum, drop=F]),
                 diag(diag(crossprod(modelCV$Model$To[, 1:OrthoLVsNum, drop=F]))^-1, ncol=OrthoLVsNum))
    }))
  loadings <- mclapply(X = 1:ntable, mc.cores = detectCores, FUN = function(i) { 
    cbind(loadings.list[[1]][[i]], loadings.list[[2]][[i]])
  })
  # Add RV coefficients in the model objects
  modelCV$RV <- RV  
  # Add the loadings in the model objects
  modelCV$Model$loadings <- loadings 
  
  # Return the final model
  return(modelCV)
}