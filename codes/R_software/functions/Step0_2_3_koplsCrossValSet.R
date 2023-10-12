#' koplsCrossValSet
#' Generates set of training/test observations for the Cross-Validation.
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
#' @param K: matrix. Kernel matrix. 
#' @param Y: matrix. Response matrix.
#' @param modelFrac: numeric. Fraction (in percent) of observations used in the 
#' training set. By default is 2/3.
#' @param type: character. Type of cross-validation wanted. It can be `nfold` 
#' for n-fold, `mccv` for Monte Carlo CV, `mccvb` for Monte Carlo class-balanced 
#' CV. Default is `nfold`.
#' @param nfold: numeric. Number of total nfold rounds (if type = 'nfold').
#' @param nfoldRound: numeric. Current nfold rounds (if type = 'nfold').
#'
#' @return
#' A list with the following entries:
#' `type`: character. Cross-validation type.
#' `nfold`: numeric. Number of nfold rounds.
#' `nfoldRound`: numeric. The current nfold round.
#' `KTrTr`: matrix. Kernel training matrix; KTrTr = <phi(Xtr),phi(Xtr)>.
#' `KTeTr`: matrix. Kernel test/training matrix; KTeTr = <phi(Xte),phi(Xtr)>.}
#' `KTeTe`: matrix. Kernel test matrix; KTeTe = <phi(Xte),phi(Xte)>.
#' `yTraining`: matrix. Y training set.
#' `yTest`: matrix. Y test set.
#' `trainingIndex`: integer. Indices of training set observations.
#' `testIndex`: numeric. Indices of test set observations.
#' `class`: character. Object class is `crossValSet`.
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
  
  # Package installation control
  if (is.element('stats', installed.packages()[,1]) == FALSE){
    install.packages('stats')
  }

  # Function loading control
  warning("Remember to load the source code for the `koplsReDummy` function.")
  
  # Define Monte-Carlos Cross Validation - class Balanced
  if(type == "mccvb"){
    # check if Y is dummy or labels
    if(base::all(base::unique(Y) == c(0, 1))){
      classVect <- koplsReDummy(Y = Y)
    } else{
      classVect <- Y
    }
    
    # Find all class labels
    uniqueClass <- unique(classVect)
    
    # For each class
    for(i in 1:length(uniqueClass)){
      # Find samples of current class
      ind <- which(classVect == uniqueClass[i])
      
      # Create a random indices
      randomSeq <- stats::rnorm(length(ind))
      rand_ind <- base::sort(randomSeq, index.return = TRUE)$ix
      ind <- ind[rand_ind]
      
      # Calculates the sample size of the training data
      trainSize <- base::floor(length(ind)*modelFrac)
      
      # Divides the sample into train and test
      trainInd <- c(trainInd, ind[1:trainSize])
      predInd <- c(predInd, ind[(trainSize+1):length(ind)])
    }
  }
  
  # Define Monte-Carlos Cross Validation
  if(type == "mccv"){
    # Create a random indices
    randomSeq <- stats::rnorm(nrow(K))
    rand_ind <- base::sort(randomSeq, index.return = TRUE)$ix
    
    # Calculates the sample size of the training data
    trainSize <- base::floor(nrow(K)*modelFrac)
    
    # Divides the sample into train and test
    trainInd <- rand_ind[1:trainSize]
    predInd <- rand_ind[(trainSize+1):nrow(K)]
  }
  
  # Define N-fold Cross Validation
  if(type == "nfold"){
    predInd <- base::seq(from = nfoldRound, to = nrow(K), by = nfold)
    trainInd <- base::setdiff(x = 1:nrow(K), y = predInd)
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
  
  # Return the final CV Set
  return(list("type" = type,
              "nfold" = nfold,
              "nfoldRound" = nfoldRound,
              "KTrTr" = KTrTr,
              "KTeTr" = KTeTr,
              "KTeTe" = KTeTe,
              "yTraining" = Y[trainInd, , drop=FALSE],
              "yTest" = Y[predInd, , drop=FALSE],
              "trainingIndex" = trainInd,
              "testIndex" = predInd,
              "class" = "crossValSet"))
}
