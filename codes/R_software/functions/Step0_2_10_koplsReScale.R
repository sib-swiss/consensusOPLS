#' koplsRescale
#' Scales a matrix based on pre-defined parameters from a scaling
#' object defined in a list (result of 'koplsScale' function).
#'
#' @param scaleS: a list containing scaling parameters 
#' @param varargin: If defined, this matrix will be scaled and returned.
#' Otherwise the original data set in the scaleS object will be scaled and 
#' returned. 
#'
#' @return
#' `scaleS`: A list containing the following entries:
#' centerType: 'mc' (mean-centering) or 'no' (no centering).
#' scaleType: 'uv' (unit variance), 'pa' (pareto) or 'no' (no scaling).
#' meanV: vector with mean values for all columns in X.
#' stdV: vector with standard deviations for all columns in X.
#' X: Scaled version of 'varargin', if defined, otherwise, scaled version of 
#' scaleS.X from input. Scaling is done according to 'centerType' and 'scaleType'.
#'
#' @examples
#'data <- matrix(c(-1.732051, 0, 1.732051, 2, 0, -2, -1.732051, 0, 1.732051), 
#'                  nrow = 3, ncol = 3)
#'scaleS <- list("centerType" = "mc", "scaleType" = "pa", "meanV" = 0, 
#'               "stdV" = 1.581139, "matrix" = data)
#'test <- koplsRescale(scaleS)
#'test
#'test$X

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
    X <- apply(X, 2, FUN = function(X){X - scaleS$meanV})
  }
  
  # Scale the matrix
  if(scaleS$scaleType == "uv"){
    X <- apply(X, 2, FUN = function(X){X/scaleS$stdV})
  }
  if(scaleS$scaleType == "pa"){
    X <- apply(X, 2, FUN = function(X){X/sqrt(scaleS$stdV)})
  }
  
  # Return the list of parameters
  return(scaleS = list("centerType" = scaleS$centerType,
                       "scaleType" = scaleS$scaleType,
                       "meanV" = mean(X), 
                       "stdV" = sd(X),
                       "X" = X))
}
