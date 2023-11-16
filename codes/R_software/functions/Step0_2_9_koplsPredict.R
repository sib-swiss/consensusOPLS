#' koplsPredict
#' Performs prediction of new samples from an existing K-OPLS model.
#' The function projects the Y-predictive and Y-orthogonal scores components 
#' to predict a value of the response matrix Y. The dimensions of the 
#' parameters is determined from the specified model.
#' 
#' # ------------------------------------------------------------------------ #
#' This file is part of the K-OPLS package, developed by Max Bylesjo, 
#' University of Umea, Judy Fonville and Mattias Rantalainen, Imperial College.
#' 
#' Copyright (c) 2007-2010 Max Bylesjo, Judy Fonville and Mattias Rantalainen 
#' 
#' This code has been extended and adapted under the terms of the GNU General 
#' Public License version 2 as published by the Free Software Foundation.
#' # ------------------------------------------------------------------------ #
#'
#' @param KteTr: matrix. The hybrid test/training kernel matrix. 
#' KteTr = <phi(Xte),phi(Xtr)>.
#' @param Ktest: matrix. The pure test kernel matrix. Ktest = <phi(Xte),phi(Xte)>.
#' @param Ktrain: matrix. The training kernel matrix (same as used in model 
#' training). Ktrain = <phi(Xtr),phi(Xtr)>.
#' @param model: list. K-OPLS model object.
#' @param nox: numeric. Number of Y-orthogonal components. If not specified, the 
#' number used during model training will be employed.
#' @param rescaleY: logical. If `TRUE`, predicted values of the
#' response (Yhat) is rescaled according to the pre-processing settings of 
#' the model. If `FALSE`, Yhat is not rescaled (default).
#'
#' @return
#' A list with the following entries:
#' `Tp`: matrix. Predicted predictive score matrix for all generations 0:'nox' 
#' of Y-orthogonal vectors.
#' `to`: vector. Predicted Y-orthogonal score vectors.
#' `T`: matrix. Predictive score matrix for the final model.
#' `KteTr`: matrix. Predictive score matrix for the final model with 'nox' 
#' Y-orthogonal vectors.
#' `EEprime`: matrix. Calculated residuals for the test kernel 'Ktest', useful 
#' e.g. for residual statistics.
#' `Yhat`: matrix. Predicted values of the response matrix.
#'
#' @examples
#' TO DO

koplsPredict <- function(KteTr, Ktest, Ktrain,
                         model, nox, rescaleY){
  # Variable format control
  if(!is.matrix(KteTr)){stop("KteTr is not a matrix.")}
  if(!is.matrix(Ktest)){stop("Ktest is not a matrix.")}
  if(!is.matrix(Ktrain)){stop("Ktrain is not a matrix.")}
  if(!is.list(model)){
    stop("model is not a list containing model parameters.")
  } else{if(model$class != "kopls"){stop("Model must be of type `kopls`.")}}
  if(!is.null(nox)){
    if(!is.numeric(nox)){stop("nox is not numeric.")}
    if(nox > model$nox){
      warning("Number of Y-orthogonal components to use is higher than in model.
              Setting number of Yorth to max in model.")
      nox <- model$nox
    }
  } else{
    stop('Number of Y-orthogonal components to use is missing.')
  }
  if(is.null(rescaleY)){
    rescaleY <- 0
  } else{if(!is.logical(rescaleY)){stop("rescaleY is not logical.")}}
  
  # Function loading control
  if (!exists("koplsCenterKTeTe", mode = "function")) {
    warning("Remember to load the source code for the `koplsCenterKTeTe` function.")
  }
  if (!exists("koplsCenterKTeTr", mode = "function")) {
    warning("Remember to load the source code for the `koplsCenterKTeTr` function.")
  }
  if (!exists("koplsRescale", mode = "function")) {
    warning("Remember to load the source code for the `koplsRescale` function.")
  }
  
  # Step1: mean centering of K matrices
  # the order of the code below is important
  KteTeMc <- Ktest
  if (model$preProc$K == "mc") {
    KteTeMc <- koplsCenterKTeTe(KteTe = Ktest, KteTr = KteTr, KtrTr = Ktrain)
  }
  KteTe <- base::matrix(data = list(NULL), nrow = model$nox+1, ncol = model$nox+1)
  KteTe[1,1][[1]] <- KteTeMc
  
  KteTrMc <- KteTr
  if (model$preProc$K == "mc") {
    KteTrMc <- koplsCenterKTeTr(KteTr, Ktrain)
  }
  KteTrTmp <- base::matrix(data = list(NULL), nrow = model$nox+1, ncol = model$nox+1)
  KteTrTmp[1,1][[1]] <- KteTrMc
  KteTr < -KteTrTmp
  
  # Initialize variables
  to <- list() ; Tp <- list()
  
  # Step2: KOPLS prediction
  ## Step2.1: for each Y-orth component
  if(nox > 0){
    for(i in 1:nox){
      ## Step2.2: Predicted predictive score matrix
      Tp[[i]] <- KteTr[i,1][[1]] %*% model$Up %*% model$Sps
      
      # Step2.3: Predicted Y-orthogonal score vectors
      to[[i]] <- (KteTr[i,i][[1]] - Tp[[i]]%*%t(model$Tp[[i]])) %*% 
        model$Tp[[i]] %*% model$co[[i]] %*% model$so[[i]]**(-1/2)
      
      # Step2.4: Normalize to
      to[[i]] <- to[[i]]/model$toNorm[[i]]
      
      # Step2.4.5: deflate KteTe (this is an EXTRA feature - not in alg. in paper)
      KteTe[i+1,i+1][[1]] <- KteTe[i,i][[1]] - 
        KteTr[i,i][[1]] %*% model$to[[i]] %*% t(to[[i]]) - 
        to[[i]] %*% t(model$to[[i]]) %*% t(KteTr[i,i][[1]]) + 
        to[[i]] %*% t(model$to[[i]]) %*% model$K[i,i][[1]] %*% model$to[[i]] %*% t(to[[i]])
      
      # Step2.5: Update KTeTr
      KteTr[i+1,1][[1]] <- KteTr[i,1][[1]] - to[[i]] %*% t(model$to[[i]]) %*% t(model$K[1,i][[1]])
      
      # Step2.6: Update KTeTr
      KteTr[i+1,i+1][[1]] <- KteTr[i,i][[1]] - 
        KteTr[i,i][[1]] %*% model$to[[i]] %*% t(model$to[[i]]) - 
        to[[i]] %*% t(model$to[[i]]) %*% model$K[i,i][[1]] + 
        to[[i]] %*% t(model$to[[i]]) %*% model$K[i,i][[1]] %*% model$to[[i]] %*% t(model$to[{i}])
    } # Step2.7: end loop
  }
  
  if (nox == 0) {
    i <- 0
  }
  
  Tp[[i+1]] <- KteTr[i+1,1][[1]] %*% model$Up %*% model$Sps
  Yhat <- Tp[[i+1]] %*% model$Bt[[i+1]] %*% t(model$Cp)
  
  if(!is.null(rescaleY)){
    if(model$preProc$Y == "no"){
      YhatRescaled <- koplsRescale(scaleS = model$preProc$paramsY, varargin = Yhat)
      Yhat <- YhatRescaled$X
    }else{
      warning("Attempted re-scale of Yhat although no pre-processing 
              parameters have been set.")
    }
  }
  
  #---- Extra stuff ----------------------------------
  #this appears to be correct - but does not match previous code...
  EEprime <- KteTe[i+1,i+1][[1]] - Tp[[i+1]] %*% t(Tp[[i+1]])
  #--------------------------------------------------
  
  # Return the list of prediction parameters
  return(list("Tp" = Tp,
              "to" = to,
              "T" = Tp[[nox+1]],
              "KteTr" = KteTr,
              "EEprime" = EEprime,
              "Yhat" = Yhat))
}