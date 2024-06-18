#' @title koplsMaxClassify
#' @description The function assign a sample to the class with maximum
#' predicted response.
#' @param X A numeric matrix of predicted responses, with class in columns and
#' sample in rows.
#'
#' @returns A vector of classes for all samples.
#'
#' @examples
#' data <- as.matrix(data.frame(x1 = c(2, 1, 5, 1),
#'                              x2 = c(7, 1, 1, 5),
#'                              x3 = c(9, 5, 4, 9),
#'                              x4 = c(3, 4, 1, 2)))
#' test <- ConsensusOPLS:::koplsMaxClassify(X = data)
#' test
#' 
#' @keywords internal
#' @noRd
#' 
koplsMaxClassify <- function(X) {
    # Variable format control
    if (!is.matrix(X)) stop("X is not a matrix.")
    
    # Search position of max value
    predClass <- colnames(X)[apply(X = X, MARGIN = 1,
                                   FUN = function(row) which.max(row))]
    
    return (predClass)
}
