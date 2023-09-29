#' koplsReDummy
#' Reconstructs a (integer) class vector from a binary (dummy) matrix.
#'
#' @param Y: a dummy matrix.
#'
#' @return 
#' `classVect`: The reconstructed integer class vector.
#'
#' @examples
#' matrix_test <- c(5, 1, 2, 3, 4, 3, 2, 4, 3, 1, 3)
#' Y <- koplsDummy(class = matrix_test, numClasses = NA)
#' Y$matrix

koplsReDummy <- function(Y){
  # Create an empty vector
  classVect <- base::rep(x = NA, times = nrow(Y))
  
  # Rebuild the vector
  for(i in 1:ncol(Y)){
    classVect[Y[, i] == 1] <- i
  }
  
  # Return a (integer) class vector from a binary matrix
  return(classVect)
}