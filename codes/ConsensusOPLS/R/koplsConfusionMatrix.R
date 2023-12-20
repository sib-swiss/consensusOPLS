#' @title koplsConfusionMatrix
#' @description Calculates a confusion matrix from classification results.
#' 
#' @param true_class vector. Indicates the true class belonging.
#' @param pred matrix. predicted class assignment. 
#'
#' @return Confusion matrix.
#'
#' @examples
#' true_class <- c(1, 2, 1, 2, 3, 3)
#' pred <- matrix(data = c(1, 2, 3, 2, 3, 3), nrow = length(true_class), ncol = 1)
#' test <- ConsensusOPLS:::koplsConfusionMatrix(true_class = true_class, 
#'                                              pred = pred)
#' test


koplsConfusionMatrix <- function(true_class, pred){
    # Variable format control
    if (!is.numeric(true_class)) stop("true_class is not numeric.")
    if (!is.matrix(pred)) stop("pred is not a matrix.")
    
    # Extract all classes
    uniqueClass <- unique(x = true_class)
    
    # Check uniqueClass format
    if (!is.numeric(uniqueClass)){
        a1 <- koplsDummy(X = true_class, numClasses = NA)
        colnames(a1) <- 1:nrow(a1)
        
        true2 <- matrix(data = NA, nrow = nrow(a1), ncol = 1)
        true_class <- apply(X = a1, MARGIN = 2,
                            FUN = function(X){
                                true2[X > 0, 1] <- colnames(X)
                            })
    }
    
    # Initialize parameters
    A <- matrix(data = 0, nrow = length(uniqueClass),
                ncol = length(uniqueClass))
    
    # For each class, find the true class indices
    indTrue <- lapply(X = uniqueClass, 
                      FUN = function(cls) which(true_class == cls))
    
    # For each class, find the corresponding prediction
    predIndex <- match(x = pred, table = uniqueClass)
    
    # Calculating the occurrences of each unique class
    pred_freqs <- sapply(X = indTrue,
                         FUN = function(indices){
                             if(length(indices) > 0){
                                 pred_classes <- predIndex[indices]
                                 pred_counts <- table(pred_classes)
                                 pred_counts/ length(indices)
                             } else{
                                 rep(x = 0, 
                                     times = length(uniqueClass))
                             }
                         })
    
    #bug here..
    A[] <- t(pred_freqs)
    
    # Return the confusion matrix
    return(A)
}
