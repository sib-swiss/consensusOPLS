#' koplsScale
#' Function for mean-centering and scaling of a matrix.
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
#' @param X: matrix. The matrix to be centered and scaled.
#' @param centerType: character. Indicates the centering method of the X matrix: 
#' mean centering (`mc`) or no centering (`no`). Default is `no` centering.
#' @param scaleType: character. Indicates the scaling method of the X matrix: 
#' `uv` for unit variance scaling (`pa`) for Pareto scaling or `no` for no 
#' scaling. Default is 'no' scaling.
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
#' Y$matrix

koplsScale <- function(X, centerType = "no", scaleType = "no"){
  # Variable format control
  if(!is.matrix(X)){stop("X is not a matrix.")}
  if(!is.character(centerType)){stop("centerType is not a character.")
  } else{ 
    if(!(centerType %in% c("mc", "no"))){
      stop("centerType must be `mc` or `no`.")}
  }
  if(!is.character(scaleType)){stop("scaleType is not a character.")
  } else{ 
    if(!(scaleType %in% c("uv", "pa", "no"))){
      stop("scaleType must be `uv`, `pa` or `no`.")}
  }
  
  # Calculation of dispersion parameters before center and scale matrix
  meanV <- base::apply(X = X, 2, FUN = function(X){mean(X)})
  stdV <- base::apply(X = X, 2, FUN = function(X){sd(X)})
  
  # Center the matrix
  if(centerType == "mc"){
    X <- base::apply(X = X, MARGIN = 2, FUN = function(X){X - mean(X)})
  }
  
  # Scale the matrix
  if(scaleType == "uv"){
    X <- base::apply(X = X, MARGIN = 2, FUN = function(X){X/sd(X)})
  }
  if(scaleType == "pa"){
    X <- base::apply(X = X, MARGIN = 2, FUN = function(X){X/sqrt(sd(X))})
  }
  
  # Return a list with all parameters
  return(list("centerType" = centerType,
              "scaleType" = scaleType,
              "meanV" = meanV,
              "stdV" = stdV,
              "matrix" = X))
}
