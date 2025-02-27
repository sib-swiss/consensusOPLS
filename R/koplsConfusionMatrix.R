#' @title koplsConfusionMatrix
#' @description Calculates a confusion matrix from classification results.
#' 
#' @param true_class vector. Indicates the true class belonging.
#' @param pred_class matrix. predicted class assignment. 
#'
#' @returns Confusion matrix.
#'
#' @examples
#' trueClass <- c(1, 2, 1, 2, 3, 3)
#' predClass <- c(1, 2, 3, 2, 3, 3)
#' ConsensusOPLS:::koplsConfusionMatrix(trueClass = trueClass, 
#'                                      predClass = predClass)
#' @keywords internal
#' @noRd
#' 
koplsConfusionMatrix <- function(trueClass, predClass) {
    return (table(trueClass, predClass))
}
