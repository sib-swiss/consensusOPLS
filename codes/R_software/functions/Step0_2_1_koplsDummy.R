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
    }
  }
  if(is.na(numClasses)){
    sample_labels <- base::length(base::unique(class))
    ncol <- sample_labels
    default <- TRUE
  } else{
    if(!is.numeric(numClasses)){stop("numClasses must be numeric.")}
    sample_labels <- base::cut(x = class, breaks = numClasses,
                               ordered_result = TRUE)
    ncol <- base::length(levels(sample_labels))
    default <- FALSE
  }
  
  # Define parameters
  labels <- base::sort(base::unique(class), decreasing = FALSE)
  
  # Matrix initialization
  dummy <- base::matrix(data = 0, nrow = base::length(class), ncol = ncol)
  
  # Search for class membership
  for(i in 1:ncol){
    if(default){
      dummy[class == labels[i], i] <- 1
      
      # And change dummy's column names
      colnames(dummy) <- labels
    } else{
      labels_cut <- cbind("lower" = as.numeric(base::sub("\\((.+),.*", "\\1", 
                                                         levels(sample_labels)[i],
                                                         fixed = FALSE)),
                          "upper" = as.numeric(base::sub("[^,]*,([^]]*)\\]", "\\1", 
                                                         levels(sample_labels)[i],
                                                         fixed = FALSE)))
      dummy[(class > labels_cut[1, "lower"] & 
               class <= labels_cut[1, "upper"]), i] <- 1
      
      # And change dummy's column names
      colnames(dummy) <- c(labels_cut[1, "lower"], labels_cut[1, "upper"])
    }
  }
  
  # Change dummy row names
  rownames(dummy) <- rownames(class)
  
  # Return a list with 2 elements
  return(list("matrix" = dummy,
              "labels_sorted" = labels))
}
