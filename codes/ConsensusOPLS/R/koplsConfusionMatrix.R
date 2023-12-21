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
#' pred_class <- c(1, 2, 3, 2, 3, 3)
#' ConsensusOPLS:::koplsConfusionMatrix(true_class = true_class, 
#'                                              pred_class = pred_class)
#' @keywords internal
koplsConfusionMatrix <- function(true_class, pred_class) {
    return (table(true_class, pred_class))
}
