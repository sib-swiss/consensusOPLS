#' koplsReDummy
#' Convert a dummy matrix into an integer vector
#'
#' @param Y: a dummy matrix 
#'
#' @return classVect
#' A (integer) class vector from a binary matrix
#'
#' @examples
#' matrix_test <- c(5, 1, 2, 3, 4, 3, 2, 4, 3, 1, 3)
#' Y <- koplsDummy(class = matrix_test, numClasses = NA)
#' Y$matrix

koplsReDummy <- function(Y){
  # Extract size of matrix
  n <- nrow(Y)
  m <- ncol(Y)
  
  # Create an empty vector
  classVect <- base::rep(x = NA, times = n)
  
  # Rebuild the vector
  for(i in 1:m){
    classVect[Y[, i] == 1] <- i
  }
  
  # Return a (integer) class vector from a binary matrix
  return(classVect)
}