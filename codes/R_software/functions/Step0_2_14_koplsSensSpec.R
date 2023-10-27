#' koplsSensSpec
#' Calculates sensitivity and specificity in a class-wise fashion.
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
#' @param trueClass: matrix. Row vector of true class assignments (template). 
#' @param predClass: matrix. Matrix (or row vector) of class assignments to be 
#' compared.
#'
#' @return
#' a list containing:
#' `classResults `: list. The results for each class.
#' `totalResults `: list. The results over all classes and averages.
#'
#' @examples
#' TO DO

koplsSensSpec <- function(trueClass, predClass){
  # Variable format control
  if(!is.data.frame(trueClass) & !is.matrix(trueClass)){
    warning("trueClass is neither a matrix nor a data.frame,
            so it was converted into a matrix.")
    trueClass <- as.matrix(trueClass)
  }
  if(!is.data.frame(predClass) & !is.matrix(predClass)){
    warning("predClass is neither a matrix nor a data.frame, 
            so it was converted into a matrix.")
    predClass <- as.matrix(predClass)
  }
  
  # Check dimensions
  if(is.vector(trueClass) & is.vector(predClass)){
    if(length(predClass) != length(trueClass)){
      stop("Template vector (trueClass) differs in length from the vector 
           (predClass) to be compared.")
    }
    
    # Forces matrix conversion for the code below
    predClass <- as.matrix(predClass)
  } else{
    if(is.vector(trueClass) & is.matrix(predClass)){
      if(nrow(predClass) != base::length(trueClass)){
        stop("Template vector (trueClass) differs in length from the matrix 
             (predClass) to be compared.")
      }
    }
  }
  
  # Function loading control
  if(!exists("koplsDummy", mode = "function")){
    warning("Remember to load the source code for the `koplsDummy` function.")
  }
  
  # To make sure dummy is correct dimension for prediction data set
  tmp1 <- base:: unique(as.vector(trueClass))
  if(ncol(trueClass) == 1){
    trueClassDummy <- koplsDummy(class = trueClass, numClasses = NA)
    ########## BUG ICI ##########
    predClassDummy <- koplsDummy(class = predClass, 
                                 numClasses = ncol(trueClassDummy$matrix))
  }
  if(length(tmp1) == 2){
    if(all(tmp1 == c(0,1))){
      trueClassDummy <- trueClass
      predClassDummy <- predClass
    }
  }
  if(nrow(predClassDummy$matrix) != nrow(trueClassDummy$matrix)){
    stop("Different number of observations in predClass and trueClass")
  }
  
  # Define classes to compare
  nclasses <- ncol(trueClassDummy)
  
  # Define TruePositifs
  TP <- base::rowSums(x = trueClass & predClass)
  # Define TrueNegatifs
  TN <- base::rowSums(x = !trueClass & !predClass)
  # Define FalsePositifs
  FN <- base::rowSums(x = trueClass & !predClass)
  # Define FalseNegatifs
  FP <- base::rowSums(x = !trueClass & predClass)
  
  # Sensitivity for each class
  sens <- TP / (TP + FN)
  # Specificity for each class
  spec <- TN / (TN + FP)
  
  # Correct missing values
  sens[is.na(sens)] <- 0
  spec[is.na(spec)] <- 0
  
  # Results for each class
  results <- list()
  for (i in 1:nclasses) {
    results[[i]] <- list(
      TP = TP[i],
      TN = TN[i],
      FP = FP[i],
      FN = FN[i],
      sens = sens[i],
      spec = spec[i],
      class = i
    )
  }
  
  # For all classes
  # Total sensibility
  sensTot <- base::sum(TP) / base::sum(TP + FN)
  sensTot[is.na(sensTot)] <- 0
  # Total specificity
  specTot <- base::sum(TN) / base::sum(TN + FP)
  specTot[is.na(specTot)] <- 0
  
  # Mean of total items
  meanSens <- base::mean(sens)
  meanSens[is.na(meanSens)] <- 0
  meanSpec <- base::mean(spec)
  meanSpec[is.na(meanSpec)] <- 0
  
  # Result for all classes
  resultsTot <- list(
    TPtot = base::sum(TP),
    FPtot = base::sum(FP),
    TNtot = base::sum(TN),
    FNtot = base::sum(FN),
    Ntot = base::sum(TP + FN),
    sensTot = sensTot,
    specTot = specTot,
    meanSens = meanSens,
    meanSpec = meanSpec
  )
  
  #Return the lists
  return(list("totalResults" = resultsTot,
              "classResults" = results))
}