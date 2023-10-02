#' koplsConfusionMatrix
#' Calculates a confusion matrix from classification results.
#'
#' @param true: Indicates the true class belonging.
#' @param pred: predicted class assignment. 
#'
#' @return
#' `A`: Confusion matrix.
#'
#' @examples
#' TO DO

koplsConfusionMatrix <- function(true, pred){
  # Variable format control
  #if(!is.vector(class)){stop("class is not a vector.")}
  if(!is.matrix(pred)){stop("pred is not a matrix.")}
  
  # Extract all classes
  uniqueClass <- unique(true)
  
  # Initialize parameters
  A <- base::matrix(0, nrow = length(uniqueClass),
                    ncol = length(uniqueClass))
  
  # For each class
  for(i in 1:length(uniqueClass)){
    # Find the indices where uniqueClass equals the current class
    indTrue <- which(true == uniqueClass[i])
    for(j in 1:length(indTrue)){
      # Find the corresponding prediction
      predIndex <- which(uniqueClass == pred[indTrue[j]])
      # Update the co-occurrence matrix
      A[i, predIndex] <- A[i, predIndex] + 1
    }
    # Normalize the row by dividing by the number of occurrences
    A[i, ] <- A[i, ] / length(indTrue)
  }
  
  # Return the confusion matrix
  return(A)
}