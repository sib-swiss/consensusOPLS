#' @title koplsScale
#' @description Function for mean-centering and scaling of a matrix.
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
#' @param X matrix. The matrix to be centered and scaled.
#' @param centerType character. Indicates the centering method of the \code{X} 
#' matrix: mean centering (\code{mc}) or no centering (\code{no}). Default is 
#' \code{no} centering.
#' @param scaleType character. Indicates the scaling method of the \code{X} 
#' matrix: \code{uv} for unit variance scaling (\code{pa}) for Pareto scaling or 
#' \code{no} for no scaling. Default is \code{no} scaling.
#'
#' @return A list containing the following entries:
#' \item{centerType}{ character. Indicates the centering method of the X matrix.}
#' \item{scaleType}{ character. Indicates the scaling method of the X matrix.}
#' \item{meanV}{ vector. Contains the mean values for all columns in X.}
#' \item{stdV}{ vector. Contains the standard deviations for all columns in X.}
#' \item{matrix}{ matrix. Original input matrix X, scaled according to 
#' \code{centerType} and \code{scaleType}.}
#'
#' @examples
#' X <- base::matrix(data = c(1,4,7, 8,4,0, 3,6,9), nrow = 3)
#' Y <- koplsScale(X = X, centerType = "mc", scaleType = "pa")
#' Y$matrix
#'  
#' @import stats
#' @keywords internal

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
  #scaleType <- match.arg(scaleType, c("uv", "pa", "no"))
  
  # Calculation of dispersion parameters before center and scale matrix
  meanV <- base::colMeans(X)
  stdV <- base::apply(X = X, 2, FUN = function(X){stats::sd(X)})
  
  # Center the matrix
  if(centerType == "mc"){
    X <- base::scale(x = X, center = TRUE, scale = FALSE)
  }
  
  # Scale the matrix
  if(scaleType == "uv"){
    X <- base::scale(x = X, center = FALSE, scale = TRUE)
  }
  if(scaleType == "pa"){
    X <- base::apply(X = X, MARGIN = 2, FUN = function(col){ col/sqrt(stats::sd(col))})
  }
  
  # Return a list with all parameters
  return(list("centerType" = centerType,
              "scaleType" = scaleType,
              "meanV" = meanV,
              "stdV" = stdV,
              "matrix" = X))
}
