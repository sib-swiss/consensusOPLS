#' koplsScale
#' Function for mean-centering and scaling of a matrix.
#'
#' @param X: the matrix to be centered and scaled.
#' @param centerType: character to indicate mean centering (`mc`) or no 
#' centering (`no`) the X matrix. Default is no centering.
#' @param scaleType: character to indicate the X matrix scaling type: `uv` for 
#' unit variance scaling (used by default), (`pa`) for Pareto scaling or `no`
#' for no scaling.
#'
#' @return
#' a list with:
#' `centerType` = 'mc' or 'no'.
#' `scaleType` = 'uv', 'pa' or 'no'.
#' `meanV` = vector with mean values for all columns in X.
#' `stdV` = vector with standard deviations for all columns in X.
#' `matrix` = Original input matrix X, scaled according to 'centerType' and 'scaleType'.
#'
#' @examples
#' X <- matrix(c(1,4,7, 8,4,0, 3,6,9), nrow = 3)
#' Y <- koplsScale(X, centerType = "mc", scaleType = "pa")
#' Y$matrix

koplsScale <- function(X, centerType = "no", scaleType = "uv"){
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
  
  # Center the matrix
  if(centerType == "mc"){
    X <- apply(X, 2, FUN = function(X){X - mean(X)})
  }
  
  # Scale the matrix
  if(scaleType == "uv"){
    X <- apply(X, 2, FUN = function(X){X/sd(X)})
  }
  if(scaleType == "pa"){
    X <- apply(X, 2, FUN = function(X){X/sqrt(sd(X))})
  }
  
  # Return a list with all parameters
  return(list("centerType" = centerType,
              "scaleType" = scaleType,
              "meanV" = mean(X),
              "stdV" = sd(X),
              "matrix" = X))
}
