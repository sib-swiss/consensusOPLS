#' koplsScaleApply
#' Applies scaling from external scaling objects on a matrix X.
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
#' @param model: list. An object containing scaling parameters.
#' @param X: matrix. The matrix to be scaled according to model parameters.
#'
#' @return
#' a list with:
#' `centerType`: character. Indicates the centering method of the X matrix.
#' `scaleType`: character. Indicates the scaling method of the X matrix.
#' `meanV`: vector. Contains the mean values for all columns in X.
#' `stdV`: vector. Contains the standard deviations for all columns in X.
#' `matrix`: matrix. Original input matrix X, scaled according to 
#' 'centerType' and 'scaleType'.
#'
#' @examples
#' X <- base::matrix(c(1,4,7, 8,4,0, 3,6,9), nrow = 3)
#' Y <- koplsScale(X, centerType = "mc", scaleType = "pa")
#' Z <- koplsScaleApply(model = Y, X = Y$matrix)
#' Z$matrix

koplsScaleApply <- function(model, X){
  # Variable format control
  if(!is.list(model)){stop("model is not a list with scaling parameters.")}
  if(!is.matrix(X)){stop("X is not a matrix.")}
  
  if(model$centerType == "mc"){
    X <- X-model$meanV
  }
  
  # Scale the matrix
  if(model$scaleType == "uv"){
    X <- X/model$stdV
  }
  if(model$scaleType == "pa"){
    X <- X/base::sqrt(model$stdV)
  }
  
  # Return a list with all parameters
  return(list("centerType" = model$centerType,
              "scaleType" = model$scaleType,
              "meanV" = model$meanV,
              "stdV" = model$stdV,
              "matrix" = X))
}
