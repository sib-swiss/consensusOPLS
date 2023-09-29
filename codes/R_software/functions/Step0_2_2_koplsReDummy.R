#' koplsReDummy
#' Reconstructs a (integer) class vector from a binary (dummy) matrix.
#'
#' @param Y: a dummy matrix.
#'
#' @return 
#' `classVect`: The reconstructed integer class vector.
#'
#' @examples
#' class_test <- c(5, 1, 2, 3, 4, 3, 2, 4, 3, 1, 3)
#' #numClasses needs to be NA for this example
#' Y <- koplsDummy(class = class_test, numClasses = NA)
#' X <- koplsReDummy(Y$matrix)

koplsReDummy <- function(Y){
  # Variable format control
  if(is.null(Y)){stop("Y must be contain matrix.")
  }else{
    test <- Y %in% c(0,1)
    for(i in 1:base::length(test)){
      if(test[i] == FALSE){
        stop("Y must contains only 0 and 1 values.")
      }
    }
  }
  
  # Create an empty vector
  classVect <- base::rep(x = NA, times = nrow(Y))
  
  # Rebuild the vector
  for(i in 1:ncol(Y)){
    classVect[Y[, i] == 1] <- i
  }
  
  # Return a (integer) class vector from a binary matrix
  return(classVect)
}
