#' koplsReDummy
#' Reconstructs a (integer) class vector from a binary (dummy) matrix.
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
#' @param Y: matrix. A dummy matrix to be transformed into a numeric vector.
#'
#' @return 
#' `classVect`: matrix. The reconstructed integer class vector.
#'
#' @examples
#' class_test <- base::matrix(c(5, 1, 2, 3, 4, 3, 2, 4, 3, 1, 3),
#'                            ncol = 1)
#' #numClasses needs to be NA for this example
#' Y <- koplsDummy(class = class_test, numClasses = NA)
#' X <- koplsReDummy(Y$matrix)

koplsReDummy <- function(Y){
  # Variable format control
  if(is.null(Y)){stop("Y must be contain  a dummy matrix.")
  }else{
    if(!is.matrix(Y)){stop("Y must be a matrix.")}
    test <- Y %in% c(0, 1)
    if(all(test) != TRUE){
      stop("Y must contains only 0 and 1 values.")
    }
  }
  
  # Create an empty vector
  classVect <- base::rep(x = NA, times = nrow(Y))
  
  # Rebuild the vector
  for(i in 1:ncol(Y)){
    classVect[Y[, i] == 1] <- i
  }
  
  # Return the reverted dummy matrix to the original vector of class labels
  return(classVect)
}
