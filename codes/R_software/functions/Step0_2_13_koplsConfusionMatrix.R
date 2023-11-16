#' koplsConfusionMatrix
#' Calculates a confusion matrix from classification results.
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
#' @param true_class: vector. Indicates the true class belonging.
#' @param pred: matrix. predicted class assignment. 
#'
#' @return
#' `A`: matrix. Confusion matrix.
#'
#' @examples
#' true_class <- c(1, 2, 1, 2, 3, 3)
#' pred <- matrix(c(1, 2, 3, 2, 3, 3), nrow = length(true_class), ncol = 1)
#' test <- koplsConfusionMatrix(true_class, pred)
#' test


koplsConfusionMatrix <- function(true_class, pred){
  # Variable format control
  if(!is.numeric(true_class)){stop("true_class is not numeric.")}
  if(!is.matrix(pred)){stop("pred is not a matrix.")}
  
  # Extract all classes
  uniqueClass <- unique(true_class)
  
  # Check uniqueClass format
  if (!is.numeric(uniqueClass)){
    # Function loading control
    if (!exists("koplsDummy", mode = "function")) {
      warning("Remember to load the source code for the `koplsDummy` function.")
    }
    
    a1 <- koplsDummy(class = true_class, numClasses = NA)
    colnames(a1) <- 1:nrow(a1)
    
    true2 <- base::matrix(data = NA, nrow = nrow(a1), ncol = 1)
    true_class <- base::apply(X = a1, MARGIN = 2,
                              FUN = function(X){
                                true2[X > 0, 1] <- colnames(X)
                                })
  }
  
  # Initialize parameters
  A <- base::matrix(0, nrow = length(uniqueClass),
                    ncol = length(uniqueClass))
  
  # For each class
  for(i in 1:length(uniqueClass)){
    # Find the indices where uniqueClass equals the current class
    indTrue <- base::which(true_class == uniqueClass[i])
    for(j in 1:length(indTrue)){
      # Find the corresponding prediction
      predIndex <- which(uniqueClass == pred[indTrue[j]])
      # Update the co-occurrence matrix
      A[i, predIndex] <- A[i, predIndex] + 1
    }
    # Normalize the row by dividing by the number of occurrences
    A[i, ] <- A[i, ] / length(indTrue)
  }
  
  # Return the confusion matrix
  return(A)
}
