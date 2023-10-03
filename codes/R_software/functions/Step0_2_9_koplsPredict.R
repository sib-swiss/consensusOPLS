#' koplsPredict
#' Performs prediction of new samples from an existing K-OPLS model.
#' The function projects the Y-predictive and Y-orthogonal scores components 
#' to predict a value of the response matrix Y. The dimensions of the 
#' parameters is determined from the specified model.
#'
#' @param KteTr: The hybrid test/training kernel matrix. KteTr = <phi(Xte),phi(Xtr)>.
#' @param Ktest: The pure test kernel matrix. Ktest = <phi(Xte),phi(Xte)>.
#' @param Ktrain: The training kernel matrix (same as used in model training). 
#' Ktrain = <phi(Xtr),phi(Xtr)>.
#' @param model: K-OPLS model object.
#' @param nox: Number of Y-orthogonal components. If not specified, the number 
#' used during model training will be employed.
#' @param rescaleY: Boolean parameter. If true, predicted values of the
#' response (Yhat) is rescaled according to the pre-processing settings of 
#' the model. If false, Yhat is not rescaled (default).
#'
#' @return
#' `modelp`: a list with the following entries:
#' `Tp`: Predicted predictive score matrix for all generations 0:'nox' of 
#' Y-orthogonal vectors.
#' `KteTr`: Predictive score matrix for the final model with 'nox' Y-orthogonal vectors.
#' `to`: Predicted Y-orthogonal score vectors.
#' `EEprime`: Calculated residuals for the test kernel 'Ktest', useful e.g. for 
#' residual statistics.
#' `Yhat`: Predicted values of the response matrix.
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
  } else{if(!is.numeric(rescaleY)){stop("rescaleY is not numeric.")}}
  
  # Step1: mean centering of K matrices
  # the order of the code below is important
  KteTeMc <- Ktest
  if (model$preProc$K == "mc") {
    KteTeMc <- koplsCenterKTeTe(Ktest, KteTr, Ktrain)
  }
  KteTe <- matrix(list(NULL), model$nox + 1, model$nox + 1)
  KteTe[1, 1] <- KteTeMc
  
  KteTrMc <- KteTr
  if (model$preProc$K == "mc") {
    KteTrMc <- koplsCenterKTeTr(KteTr, Ktrain)
  }
  KteTr <- matrix(list(NULL), model$nox + 1, model$nox + 1)
  KteTr[1, 1] <- KteTrMc
  
  # Initialize 'to' variable
  to <- NULL
  
  # Step2: KOPLS prediction
  ## Step2.1: for each Y-orth component
  for(i in 1:nox){
    ## Step2.2: Predicted predictive score matrix
    Tp[i] <- KteTr[i,1] %*% model$Up %*% model$Sp**(-1/2)
    
    # Step2.3: Predicted Y-orthogonal score vectors
    to[i] <- (KteTr[i,i] - Tp[i]%*%t(model$Tp[i])) %*% 
      model$Tp[i] %*% model$co[i] %*% model$so[i]**(-1/2)
    
    # Step2.4: Normalize to
    to[i] <- to[i]/model$toNorm[i]
    
    # Step2.4.5: deflate KteTe (this is an EXTRA feature - not in alg. in paper)
    KteTe[i+1,i+1] <- KteTe[i,i] - KteTr[i,i] %*% model$to[i] %*% t(to[i]) - 
      to[i] %*% t(model$to[i]) %*% t(KteTr[i,i]) + 
      to[i] %*% t(model$to[i]) %*% model$K[i,i] %*% model$to[i] %*% t(to[i])
    
    # Step2.5: Update KTeTr
    KteTr[i+1,1] <- KteTr[i,1] - to[i] %*% t(model$to[i]) %*% t(model$K[1,i])
    
    # Step2.6: Update KTeTr
    KteTr[i+1,i+1] <- KteTr[i,i] - KteTr[i,i] %*% model$to[i] %*% t(model$to[i]) - 
      to[i] %*% t(model$to[i]) %*% model$K[i,i] + 
      to[i] %*% t(model$to[i]) %*% model$K[i,i] %*% model$to[i] %*% t(model$to[i])
  } # Step2.7: end loop
  
  if (nox == 0) {
    i <- 0
  }
  
  Tp[i+1] <- KteTr[i+1,1] %*% model$Up %*% model$Sp**(-1/2)
  Yhat <- Tp[i+1] %*% model$Bt[i+1] %*% t(model$Cp)
  
  if(!is.null(rescaleY)){
    if(model$preProc$Y == "no"){
      YhatRescaled <- koplsRescale(model$preProc$paramsY, Yhat)
      Yhat <- YhatRescaled$X
    }else{
      warning('Attempted re-scale of Yhat although no pre-processing 
              parameters have been set.')
    }
  }
  
  #---- Extra stuff ----------------------------------
  #this appears to be correct - but does not match previous code...
  EEprime <- KteTe[i+1,i+1] - Tp[i+1] %*% t(Tp[i+1])
  #--------------------------------------------------
  
  # Return the list of prediction parameters
  return(list("Tp" = Tp,
              "to" = to,
              "KteTr" = KteTr,
              "EEprime" = EEprime,
              "Yhat" = Yhat))
}