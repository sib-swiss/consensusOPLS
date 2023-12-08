#' @title koplsReDummy
#' @description Reconstructs a (integer) class vector from a binary (dummy) 
#' matrix.
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
#' @param Y matrix. A dummy matrix to be transformed into a numeric vector.
#'
#' @return The reconstructed integer class vector.
#'
#' @examples
#' class <- base::matrix(c(5, 1, 2, 3, 4, 3, 2, 4, 3, 1, 3), ncol = 1)
#' #numClasses needs to be missing for this example
#' Y <- koplsDummy(class = class, numClasses = NA)
#' X <- koplsReDummy(Y)
#' X
#' 
#' @keywords internal

koplsReDummy <- function(Y){
  # Variable format control
  if(is.null(Y)){stop("Y must be contain a dummy matrix.")
  }else{
    if(!is.matrix(Y)){stop("Y must be a matrix.")}
    if(base::all(Y %in% c(0, 1)) != TRUE){
      stop("Y must contains only 0 and 1 values.")
    }
  }
  
  # Rebuild the vector
  X <- base::apply(X = Y, MARGIN = 1, 
                   FUN = function(X) base::which(X == 1))
  
  # Return the reverted dummy matrix to the original vector of class labels
  return(X)
}
