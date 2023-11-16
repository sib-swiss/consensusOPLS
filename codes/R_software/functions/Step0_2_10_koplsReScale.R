#' koplsRescale
#' Scales a matrix based on pre-defined parameters from a scaling
#' object defined in a list (result of 'koplsScale' function).
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
#' @param scaleS: list. It contains scaling parameters 
#' @param varargin: matrix. If defined, this matrix will be scaled and returned.
#' Otherwise the original data set in the scaleS object will be scaled and 
#' returned. 
#'
#' @return
#' A list containing the following entries:
#' `centerType`: character. 'mc' (mean-centering) or 'no' (no centering).
#' `scaleType`: character. 'uv' (unit variance), 'pa' (pareto) or 'no' (no scaling).
#' `meanV`: vector. It contains the mean values for all columns in X.
#' `stdV`: vector. It contains the standard deviations for all columns in X.
#' `X`: matrix. Scaled version of 'varargin', if defined, otherwise, scaled 
#' version of 'scaleS$matrix' from input. Scaling is done according to 
#' 'centerType' and 'scaleType'.
#'
#' @examples
#' data <- base::matrix(data = c(-1.732051, 0, 1.732051, 2, 0, 
#'                               -2, -1.732051, 0, 1.732051), 
#'                      nrow = 3, ncol = 3)
#' scaleS <- list("centerType" = "mc", "scaleType" = "pa", "meanV" = 0, 
#'                "stdV" = 1.581139, "matrix" = data)
#' test <- koplsRescale(scaleS)
#' test
#' test$X

koplsRescale <- function(scaleS, varargin = NULL){
  # Variable format control
  if(!is.list(scaleS)){stop("scaleS must be a list (result of 'koplsScale()').")}
  if(!is.null(varargin)){
    if(!is.matrix(varargin)){stop("varargin must be a matrix.")}
    X <- varargin
  } else{
    X <- scaleS$matrix
  }
  
  # Center the matrix
  if(scaleS$centerType == "mc"){
    X <- X+scaleS$meanV
  }
  
  # Scale the matrix
  if(scaleS$scaleType == "uv"){
    X <- X*scaleS$stdV
  }
  if(scaleS$scaleType == "pa"){
    X <- X*sqrt(scaleS$stdV)
  }
  
  # Return the list of parameters
  return(scaleS = list("centerType" = "no",
                       "scaleType" = "no",
                       "meanV" = scaleS$meanV, 
                       "stdV" = scaleS$stdV,
                       "X" = X))
}
