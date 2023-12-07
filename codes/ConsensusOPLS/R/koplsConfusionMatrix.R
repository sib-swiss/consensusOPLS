#' @title koplsConfusionMatrix
#' @description Calculates a confusion matrix from classification results.
#' 
#' @param true_class vector. Indicates the true class belonging.
#' @param pred matrix. predicted class assignment. 
#'
#' @return
#' \item{A}{ matrix. Confusion matrix.}
#'
#' @examples
#' true_class <- c(1, 2, 1, 2, 3, 3)
#' pred <- base::matrix(c(1, 2, 3, 2, 3, 3), nrow = length(true_class), ncol = 1)
#' test <- ConsensusOPLS:::koplsConfusionMatrix(true_class, pred)
#' test


koplsConfusionMatrix <- function(true_class, pred){
  # Variable format control
  if(!is.numeric(true_class)){stop("true_class is not numeric.")}
  if(!is.matrix(pred)){stop("pred is not a matrix.")}
  
  # Extract all classes
  uniqueClass_true <- base::unique(x = true_class)
  
  # Check uniqueClass format
  if (!is.numeric(uniqueClass_true)){
    a1 <- ConsensusOPLS:::koplsDummy(X = true_class, numClasses = NA)
    colnames(a1) <- 1:nrow(a1)
    
    true2 <- base::matrix(data = NA, nrow = nrow(a1), ncol = 1)
    true_class <- base::apply(X = a1, MARGIN = 2,
                              FUN = function(X){
                                true2[X > 0, 1] <- colnames(X)
                                })
  }
  
  # Initialize parameters
  A_true <- base::matrix(data = 0, nrow = length(uniqueClass_true),
                    ncol = length(uniqueClass_true))
  
  # For each class
  for(i in 1: base::length(uniqueClass_true)){
    # Find the indices where uniqueClass equals the current class
    indTrue_true <- base::which(true_class == uniqueClass_true[i])
    for(j in 1: base::length(indTrue_true)){
      # Find the corresponding prediction
      predIndex_true <- base::which(uniqueClass_true == pred[indTrue_true[j]])
      # Update the co-occurrence matrix
      A_true[i, predIndex_true] <- A_true[i, predIndex_true] + 1
    }
    # Normalize the row by dividing by the number of occurrences
    A_true[i, ] <- A_true[i, ] / base::length(indTrue_true)
  }
  
  # Return the confusion matrix
  return(A)
}
