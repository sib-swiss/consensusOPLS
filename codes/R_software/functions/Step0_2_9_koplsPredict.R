#' @title koplsPredict
#' @description Performs prediction of new samples from an existing K-OPLS model.
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
#' @param KteTr matrix. The hybrid test/training kernel matrix. 
#' KteTr = <phi(Xte),phi(Xtr)>.
#' @param Ktest matrix. The pure test kernel matrix. Ktest = <phi(Xte),phi(Xte)>.
#' @param Ktrain matrix. The training kernel matrix (same as used in model 
#' training). Ktrain = <phi(Xtr),phi(Xtr)>.
#' @param model list. K-OPLS model object.
#' @param nox numeric. Number of Y-orthogonal components. If not specified, the 
#' number used during model training will be employed.
#' @param rescaleY logical. If \code{TRUE}, predicted values of the
#' response (Yhat) is rescaled according to the pre-processing settings of 
#' the model. If \code{FALSE}, Yhat is not rescaled (default).
#'
#' @return A list with the following entries:
#' \item{Tp}{ matrix. Predicted predictive score matrix for all generations 0: 
#' \code{nox} of Y-orthogonal vectors.}
#' \item{to}{ vector. Predicted Y-orthogonal score vectors.}
#' \item{T}{ matrix. Predictive score matrix for the final model.}
#' \item{KteTr}{ matrix. Predictive score matrix for the final model with 
#' \code{nox} Y-orthogonal vectors.}
#' \item{EEprime}{ matrix. Calculated residuals for the test kernel \code{Ktest}, 
#' useful e.g. for residual statistics.}
#' \item{Yhat}{ matrix. Predicted values of the response matrix.}
#'
#' @examples
#' Xte <- base::matrix(data = stats::rnorm(n = 20), ncol=5)
#' Xtr <- base::matrix(data = stats::rnorm(n = 25), ncol=5)
#' KteTe <- koplsKernel(X1 = Xte, X2 = Xte, 
#'                      Ktype='g', params=c(sigma=1.0))
#' KteTr <- koplsKernel(X1 = Xte, X2 = Xtr, 
#'                      Ktype='g', params=c(sigma=1.0))
#' KtrTr <- koplsKernel(X1 = Xtr, X2 = Xtr, 
#'                      Ktype='g', params=c(sigma=1.0))
#' 
#' Y <- base::matrix(data = stats::rnorm(n = 15), nrow = 5)
#' A <- 2
#' nox <- 4
#' preProcK <- "mc"
#' preProcY <- "mc"
#' model <- koplsModel(K = KtrTr, Y = Y, A = A, nox = nox, 
#'                     preProcK = preProcK, preProcY = preProcY)
#' pred <- koplsPredict(KteTr = KteTr, Ktest = KteTe, 
#'                      Ktrain = KtrTr, model = model, nox = nox,
#'                      rescaleY = FALSE)
#' pred
#' 
#' @keywords internal

koplsPredict <- function(KteTr, Ktest, Ktrain,
                         model, nox, rescaleY){
  # Variable format control
  if (!is.matrix(KteTr) || !is.matrix(Ktest) || !is.matrix(Ktrain)) {
    stop("One or more kernel inputs are not matrices.")
  }
  if(!is.list(model)){
    stop("model is not a list containing model parameters.")
  } else{if(model$Unique_params$class != "kopls"){stop("Model must be of type `kopls`.")}}
  if(!is.null(nox)){
    if(!is.numeric(nox)){stop("nox is not numeric.")}
    if(nox > model$Unique_params$nox){
      warning("Number of Y-orthogonal components to use is higher than in model.
              Setting number of Yorth to max in model.")
      nox <- model$Unique_params$nox
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
  if (model$Unique_params$preProcK == "mc") {
    KteTeMc <- koplsCenterKTeTe(KteTe = Ktest, KteTr = KteTr, KtrTr = Ktrain)
  }
  KteTe <- base::matrix(data = list(NULL), nrow = model$Unique_params$nox+1, 
                        ncol = model$Unique_params$nox+1)
  KteTe[1,1][[1]] <- KteTeMc
  
  KteTrMc <- KteTr
  if (model$Unique_params$preProcK == "mc") {
    KteTrMc <- koplsCenterKTeTr(KteTr = KteTr, KtrTr = Ktrain)
  }
  KteTrTmp <- base::matrix(data = list(NULL), nrow = model$Unique_params$nox+1, 
                           ncol = model$Unique_params$nox+1)
  KteTrTmp[1,1][[1]] <- KteTrMc
  KteTr <- KteTrTmp
  
  # Initialize variables
  to <- list() ; Tp <- list()
  
  # Step2: KOPLS prediction
  ## Step2.1: for each Y-orth component
  if(nox > 0){
    for(i in 1:nox){
      ## Step2.2: Predicted predictive score matrix
      Tp[[i]] <- base::crossprod(x = t(KteTr[i,1][[1]]),
                                 y = base::tcrossprod(x = model$Up, 
                                                      y = t(model$Sps)))
      
      # Step2.3: Predicted Y-orthogonal score vectors
      to[[i]] <- base::crossprod(x = t((KteTr[i,i][[1]] - 
                                          base::tcrossprod(x = Tp[[i]], 
                                                           y = model$Tp[[i]]))),
                                 y = base::tcrossprod(x = model$Tp[[i]], 
                                                      y = base::tcrossprod(x = t(base::sqrt(model$so[[i]])), 
                                                                           y = t(model$co[[i]]))))
      
      # Step2.4: Normalize to
      to[[i]] <- to[[i]]/model$toNorm[[i]]
      
      # Step2.4.5: deflate KteTe (this is an EXTRA feature - not in alg. in paper)
      KteTe[i+1,i+1][[1]] <- KteTe[i,i][[1]] - 
        base::crossprod(x = t(KteTr[i,i][[1]]), 
                        y = base::tcrossprod(x = model$to[[i]], 
                                             y = to[[i]])) - 
        base::crossprod(x = t(to[[i]]), 
                        y = base::crossprod(x = model$to[[i]], 
                                            y = t(KteTr[i,i][[1]]))) +
        base::crossprod(x = t(to[[i]]),
                        y = base::crossprod(x = model$to[[i]],
                                            y = base::crossprod(x = t(model$K[i,i][[1]]),
                                                                y = base::tcrossprod(x = model$to[[i]],
                                                                                     y = to[[i]]))))
      
      # Step2.5: Update KTeTr
      KteTr[i+1,1][[1]] <- KteTr[i,1][[1]] - 
        base::crossprod(x = t(to[[i]]),
                        y = base::crossprod(x = model$to[[i]], 
                                            y = t(model$K[1,i][[1]])))
      
      # Step2.6: Update KTeTr
      KteTr[i+1,i+1][[1]] <- KteTr[i,i][[1]] - 
        base::crossprod(x = t(KteTr[i,i][[1]]),
                        y = base::tcrossprod(model$to[[i]])) - 
        base::crossprod(x = t(to[[i]]),
                        y = base::crossprod(x = model$to[[i]], 
                                            y = model$K[i,i][[1]])) + 
        base::crossprod(x = t(to[[i]]), 
                        y = base::crossprod(x = model$to[[i]], 
                                            y = base::crossprod(x = t(model$K[i,i][[1]]), 
                                                                y = base::tcrossprod(model$to[[i]]))))
    } # Step2.7: end loop
  }
  
  if (nox == 0) {
    i <- 0
  }
  
  Tp[[i+1]] <- base::crossprod(x = t(KteTr[i+1,1][[1]]),
                               y = base::tcrossprod(model$Up, t(model$Sps)))
  Yhat <- base::crossprod(x = t(Tp[[i+1]]),
                          y = base::tcrossprod(model$Bt[[i+1]], model$Cp))
  
  if(!is.null(rescaleY)){
    if(model$Unique_params$preProcY == "no"){
      if(base::length(model$preProc$paramsY) == 1 & model$preProc$paramsY == "no"){
        scaleParams <- list()
        scaleParams$centerType <- "no"
        scaleParams$scaleType <- "no"
      } else{
        scaleParams <- model$preProc$paramsY
      }
      YhatRescaled <- koplsRescale(scaleS = scaleParams, varargin = Yhat)
      Yhat <- YhatRescaled$X
    } else{
      warning("Attempted re-scale of Yhat although no pre-processing 
              parameters have been set.")
    }
  }
  
  #---- Extra stuff ----------------------------------
  EEprime <- KteTe[i+1,i+1][[1]] - base::tcrossprod(Tp[[i+1]])
  #--------------------------------------------------
  
  # Return the list of prediction parameters
  return(list("Tp" = Tp,
              "to" = to,
              "T" = Tp[[nox+1]],
              "KteTr" = KteTr,
              "EEprime" = EEprime,
              "Yhat" = Yhat))
}
