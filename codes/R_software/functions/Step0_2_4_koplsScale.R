koplsScale <- function(X, centerType, scaleType){
  # Mean-centering and scaling matrix
  # X the matrix to be centered and scaled
  # centerType for mean (or no) centering
  # scaleType for unit variance scaling, Pareto scaling of no scaling
  
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
  
  # Define variables
  m <- nrow(X)
  
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
  
  # Return 
  return(list("centerType" = centerType,
              "scaleType" = scaleType,
              "meanV" = mean(X),
              "stdV" = sd(X),
              "matrix" = X))
  # example:
  # X <- matrix(c(1,4,7, 8,4,0, 3,6,9), nrow=3)
}