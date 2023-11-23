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
#' @param numClasses numeric. Pre-defined the number of classes in the output.
#' By default (if numClasses is missing), is equal to the length of \code{class}
#' vector. Otherwise, the class vector is split into \code{numclasse} pieces.

#'
#' @return
#' \item{X}{matrix. Results of the dummy matrix, with rows corresponding to
#' observations and columns to classes. Each element in matrix is either one
#' (observation belongs to class) or zero (observation does not belong to class).}
#'
#' @examples
#' class <- matrix(c(5, 1, 2, 3, 4, 3, 2, 4, 3, 1, 3), ncol = 1)
#' Y <- ConsensusOPLS:::koplsDummy(X = class, numClasses = NA)
#' Y
#'
#' @keywords internal

koplsDummy <- function(X, numClasses = NA) {
  # Variable format control
  if (is.null(X)) {
    stop("X must be contain an integer vector.")
  } else{
    if (!is.matrix(X)) {
      stop("X must be a matrix.")
    } else{
      if (ncol(X) != 1) {
        stop("X must be a matrix with only 1 column.")
      }
      X <- as.vector(X)
    }
  }
  if (is.na(numClasses)) {
    labels <- sort(unique(X))
  } else{
    if (!is.numeric(numClasses)) {
      stop("numClasses must be numeric.")
    }
    labels <- 1:numClasses
  }
  
  # Search for class membership
  dummy <- t(sapply(
    X = X,
    FUN = function(class) {
      tmp <- numeric(length = length(labels))
      tmp[match(class, labels)] <- 1
      return(tmp)
    }
  ))
  
  # Change dummy's column names
  colnames(dummy) <- as.character(labels)
  
  # Return the matrix
  return(dummy)
}



#' @title koplsReDummy
#' @description Reconstructs a (integer) class vector from a binary (dummy) 
#' matrix. This function is the inverse of \code{koplsDummy} function.
#'
#' @param Y matrix. A dummy matrix to be transformed into a numeric vector.
#'
#' @return
#' \item{classVect}{ matrix. The reconstructed integer class vector.}
#'
#' @examples
#' class <- matrix(c(5, 1, 2, 3, 4, 3, 2, 4, 3, 1, 3), ncol = 1)
#' Y <- ConsensusOPLS:::koplsDummy(X = class, numClasses = NA)
#' X <- ConsensusOPLS:::koplsReDummy(Y)
#' X
#'
#' @keywords internal

koplsReDummy <- function(Y) {
  # Variable format control
  if (is.null(Y)) {
    stop("Y must be contain  a dummy matrix.")
  } else{
    if (!is.matrix(Y)) {
      stop("Y must be a matrix.")
    }
    if (all(Y %in% c(0, 1)) != TRUE) {
      stop("Y must contains only 0 and 1 values.")
    }
  }
  
  # Rebuild the vector
  classVect <- apply(
    X = Y,
    MARGIN = 1,
    FUN = function(X)
      which(X == 1)
  )
  
  # Return the reverted dummy matrix to the original vector of class labels
  return(classVect)
}
