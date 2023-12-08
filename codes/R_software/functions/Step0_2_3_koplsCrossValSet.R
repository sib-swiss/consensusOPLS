#' @title koplsCrossValSet
#' @description Generates set of training/test observations for the 
#' Cross-Validation.
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
#' @param K matrix. Kernel matrix. 
#' @param Y matrix. Response matrix.
#' @param modelFrac numeric. Fraction (in percent) of observations used in the 
#' training set. By default is 2/3.
#' @param type character. Type of cross-validation wanted. It can be \code{nfold} 
#' for n-fold, \code{mccv} for Monte Carlo CV, \code{mccvb} for Monte Carlo 
#' class-balanced CV. Default is \code{nfold}.
#' @param nfold numeric. Number of total nfold rounds (if type = \code{nfold}).
#' @param nfoldRound numeric. Current nfold rounds (if type = \code{nfold}).
#'
#' @return A list with the following entries:
#' \item{CV_param}{ data frame. It contains \code{type} a character for the 
#' Cross-validation type, \code{nfold} the total number of nfold, 
#' \code{nfoldRound} the current nfold round and \code{class} the object class
#' wich is \code{crossValSet}.}
#' \item{KTrTr}{ matrix. Kernel training matrix; KTrTr = <phi(Xtr),phi(Xtr)>.}
#' \item{KTeTr}{ matrix. Kernel test/training matrix; KTeTr = <phi(Xte),phi(Xtr)>.}
#' \item{KTeTe}{ matrix. Kernel test matrix; KTeTe = <phi(Xte),phi(Xte)>.}
#' \item{yTraining}{ matrix. Y training set.}
#' \item{yTest}{ matrix. Y test set.}
#' \item{trainingIndex}{ integer. Indices of training set observations.}
#' \item{testIndex}{ numeric. Indices of test set observations.}
#'
#' @examples
#' Y <- base::matrix(1:100, nrow = 10)
#' K <- base::matrix(1:100, nrow = 10)  
#' type <- "nfold"  
#' modelFrac <- 0.7  
#' nfold <- 5  
#' nfoldRound <- 1 
#' test <- koplsCrossValSet(K = K, Y = Y, modelFrac = modelFrac,
#'                          type = type, nfold = nfold, nfoldRound = nfoldRound) 
#' test                          
#'
#' @keywords internal                          
#' @import parallel

koplsCrossValSet <- function(K, Y, modelFrac = 2/3, type = "nfold", 
                             nfold = NA, nfoldRound = NA){
  # Variable format control
  if(!is.matrix(K)){stop("K is not a matrix.")}
  if(!is.matrix(Y)){stop("Y is not a matrix.")}
  if(!is.numeric(modelFrac)){stop("modelFrac is not numeric.")}
  if(!is.character(type)){stop("type is not a character.")
  } else{
    if(!(type %in% c("nfold", "mccv", "mccvb"))){
      stop("type must be `nfold`, `mccv` or `mccvb`.")
    }
  }
  if(!is.na(nfold)){
    if(!is.numeric(nfold)){stop("nfold is not a number.")}
  } else{
    if(type != "nfold"){
      warning("type is not nfold, nfold is not defined (missing argument).")}
  }
  if(!is.na(nfoldRound)){
    if(!is.numeric(nfoldRound)){stop("nfoldRound is not a number.")}
  } else{
    if(type != "nfold"){
      warning("type is not nfold, nfoldRound is not defined (missing argument).")
    }
  }
  
  # Define Monte-Carlos Cross Validation - class Balanced
  if(type == "mccvb"){
    # check if Y is dummy or labels
    if(base::all(base::unique(x = Y) == c(0, 1))){
      classVect <- koplsReDummy(Y = Y)
    } else{
      classVect <- Y
    }
    
    # Find all class labels
    uniqueClass <- base::unique(x = classVect)
    
    # Find samples of each class
    indList <- parallel::mclapply(X = uniqueClass, 
                                  mc.cores = detectCores(),
                                  FUN = function(i){
                                    ind <- base::which(classVect == i)
                                    rand_ind <- base::sample(x = ind)
                                    trainSize <- base::floor( base::length(ind)*modelFrac)
                                    c(rand_ind[1:trainSize], 
                                      rand_ind[(trainSize+1): base::length(ind)])
                                  })
    
    # Combine indices for all classes
    trainInd <- base::unlist(indList)
    predInd <- base::setdiff(x = base::seq_len(nrow(Y)), y = trainInd)
  }
  
  # Define Monte-Carlos Cross Validation
  if(type == "mccv"){
    # Create a random indices
    rand_ind <- base::sample( base::seq_len(nrow(K)))
    
    # Calculates the sample size of the training data
    trainSize <- base::floor(nrow(K)*modelFrac)
    
    # Divides the sample into train and test
    trainInd <- rand_ind[1:trainSize]
    predInd <- rand_ind[(trainSize+1):nrow(K)]
  }
  
  # Define N-fold Cross Validation
  if(type == "nfold"){
    predInd <- base::seq(from = nfoldRound, to = nrow(Y), by = nfold)
    trainInd <- base::setdiff(x = base::seq_len(nrow(Y)), y = predInd)
  }
  
  # Construct Kernel/Y matrices for training/test
  if(nrow(K) == ncol(K)){
    KTrTr <- K[trainInd, trainInd, drop=FALSE]
    KTeTr <- K[predInd, trainInd, drop=FALSE]
    KTeTe <- K[predInd, predInd, drop=FALSE]
  } else{
    KTrTr <- NA
    KTeTr <- NA
    KTeTe <- NA
  }
  
  # Group parameters in data.frame
  CV_param <- base::data.frame("type" = type,
                               "nfold" = nfold,
                               "nfoldRound" = nfoldRound,
                               "class" = "crossValSet")
  
  # Return the final CV Set
  return(list("CV_param" = CV_param,
              "KTrTr" = KTrTr,
              "KTeTr" = KTeTr,
              "KTeTe" = KTeTe,
              "yTraining" = Y[trainInd, , drop=FALSE],
              "yTest" = Y[predInd, , drop=FALSE],
              "trainingIndex" = trainInd,
              "testIndex" = predInd))
}
