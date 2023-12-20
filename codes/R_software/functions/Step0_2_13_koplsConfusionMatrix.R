#' @title koplsConfusionMatrix
#' @description Calculates a confusion matrix from classification results.
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
#' @param true_class vector. Indicates the true class belonging.
#' @param pred matrix. predicted class assignment. 
#'
#' @return
#' \item{A}{ matrix. Confusion matrix.}
#'
#' @examples
#' true_class <- c(1, 2, 1, 2, 3, 3)
#' pred <- base::matrix(data = c(1, 2, 3, 2, 3, 3), 
#'                      nrow = length(true_class), ncol = 1)
#' test <- koplsConfusionMatrix(true_class = true_class, pred = pred)
#' test


koplsConfusionMatrix <- function(true_class, pred){
  # Variable format control
  if(!is.numeric(true_class)){stop("true_class is not numeric.")}
  if(!is.matrix(pred)){stop("pred is not a matrix.")}
  
  # Extract all classes
  uniqueClass <- base::unique(x = true_class)
  
  # Check uniqueClass format
  if (!is.numeric(uniqueClass)){
    # Function loading control
    if (!exists("koplsDummy", mode = "function")) {
      warning("Remember to load the source code for the `koplsDummy` function.")
    }
    
    a1 <- koplsDummy(X = true_class, numClasses = NA)
    colnames(a1) <- 1:nrow(a1)
    
    true2 <- base::matrix(data = NA, nrow = nrow(a1), ncol = 1)
    true_class <- base::apply(X = a1, MARGIN = 2,
                              FUN = function(X){
                                true2[X > 0, 1] <- colnames(X)
                              })
  }
  
  # Initialize parameters
  A <- base::matrix(data = 0, nrow = base::length(uniqueClass),
                    ncol = base::length(uniqueClass))
  
  # For each class, find the true class indices
  indTrue <- base::lapply(X = uniqueClass, 
                          FUN = function(cls) which(true_class == cls))
  
  # For each class, find the corresponding prediction
  predIndex <- base::match(x = pred, table = uniqueClass)
  
  # Calculating the occurrences of each unique class
  pred_freqs <- base::sapply(X = indTrue,
                             FUN = function(indices){
                               if(length(indices) > 0){
                                 pred_classes <- predIndex[indices]
                                 pred_counts <- base::table(pred_classes)
                                 pred_counts/ base::length(indices)
                               } else{
                                 base::rep(x = 0, 
                                           times = base::length(uniqueClass))
                               }
                             })
  
  #bug here..
  A[] <- t(pred_freqs)
  
  # Return the confusion matrix
  return(A)
}
