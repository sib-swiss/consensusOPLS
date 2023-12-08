#' @title KoplsDummy
#' @description Transform a vector into a dummy (binary) matrix.
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
#' @param X numeric vector with classes to define.
#' @param numClasses numeric. Pre-defined number of classes in the output.
#' By default, numClasses, number of unique values in \code{X}. Otherwise,
#' classes 1,...,numClasses are considered.
#'
#' @return
#' \item{X}{matrix. Results of the dummy matrix, with rows corresponding to 
#' observations and columns to classes. Each element in matrix is either one 
#' (observation belongs to class) or zero (observation does not belong to class).}
#'
#' @examples
#' class <- base::matrix(c(5, 1, 2, 3, 4, 3, 2, 4, 3, 1, 3), ncol = 1)
#' Y <- koplsDummy(X = class, numClasses = NA)
#' Y
#' 
#' @keywords internal

koplsDummy <- function(X, numClasses = NA){
  # Variable format control
  if(is.null(X)){stop("X must be contain an integer vector.")
  } else{
    if(!is.matrix(X)){
      stop("X must be a matrix.")
    } else{
      if(ncol(X) != 1){
        stop("X must be a matrix with only 1 column.")
      }
      X <- as.vector(X)
    }
  }
  if(is.na(numClasses)){
    labels <- base::sort(x = base::unique(X))
  } else{
    if(!is.numeric(numClasses)){stop("numClasses must be numeric.")}
    labels <- 1:numClasses
  }
  
  # Search for class membership
  dummy <- t( base::sapply(X = X,
                           FUN = function(class){
                             tmp <- numeric(length = base::length(labels))
                             tmp[base::match(class, labels)] <- 1
                             return(tmp)
                           }))
  
  # Change dummy's column names
  colnames(dummy) <- as.character(labels)
  
  # Return the matrix
  return(dummy)
}
