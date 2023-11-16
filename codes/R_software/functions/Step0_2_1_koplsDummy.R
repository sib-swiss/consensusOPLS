#' KoplsDummy
#' Transform a vector into a dummy (binary) matrix.
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
#' @param class: matrix with 1 column. Integer vector with classes to define.
#' @param numClasses: numeric. Pre-defined the number of classes in the output. 
#' By default (if numClasses is missing), is equal to the length of `class` 
#' vector. Otherwise, the class vector is split into `numclasses` pieces.
#'
#' @return
#' `matrix`: matrix. Results of the dummy matrix, with rows corresponding to 
#' observations and columns to classes. Each element in matrix is either one 
#' (observation belongs to class) or zero (observation does not belong to class).
#' `labels_sorted`: numerical vector. The class labels that are found in class 
#' in sorted order.
#'
#' @examples
#' class <- base::matrix(c(5, 1, 2, 3, 4, 3, 2, 4, 3, 1, 3),
#'                       ncol = 1)
#' numClasses <- 2
#' Y <- koplsDummy(class = class, numClasses = numClasses)
#' Y$matrix
#' Y$labels_sorted

koplsDummy <- function(class, numClasses = NA){
  # Variable format control
  if(is.null(class)){stop("class must be contain an integer vector.")
  } else{
    if(!is.matrix(class)){
      stop("class must be a matrix.")
    } else{
      if(ncol(class) != 1){
        stop("class must be a matrix with only 1 column.")
      }
      class <- as.vector(class)
    }
  }
  if(is.na(numClasses)){
    labels <- base::sort(x = base::unique(class))
  } else{
    if(!is.numeric(numClasses)){stop("numClasses must be numeric.")}
    labels <- 1:numClasses
  }
  
  # Matrix initialization
  dummy <- base::matrix(data = 0, nrow = base::length(class), 
                        ncol = base::length(labels))
  
  
  for(i in 1: base::length(labels)){
    # Search for class membership
    dummy[which(class == labels[i]), i] <- 1
  }
  
  # Change dummy's row names
  rownames(dummy) <- rownames(class)
  # Change dummy's column names
  colnames(dummy) <- labels
  
  # Return a list with 2 elements
  return(list("matrix" = dummy,
              "labels_sorted" = labels))
}
