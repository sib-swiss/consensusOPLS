#' @title koplsBasicClassify
#' @description Classification function that assesses class belonging of a 
#' predicted response in \code{data} based on a fixed threshold \code{k}.
#'
#' @param X matrix. It contains the predicted response matrix Y,
#' where columns denote classes and rows observations.
#' @param k numeric. Threshold value used to assign class categories.
#'
#' @return Predicted classification.
#'
#' @examples
#' data <- as.matrix(data.frame(x1 = c(2, 1, 5, 1),
#'                              x2 = c(7, 1, 1, 5),
#'                              x3 = c(9, 5, 4, 9),
#'                              x4 = c(3, 4, 1, 2)))
#' k <- 4
#' test <- ConsensusOPLS:::koplsBasicClassify(X = data, k = k)
#' test
#' 
#' @keywords internal
#' 
koplsBasicClassify <- function(X, k) {
    # Variable format control
    if (!is.matrix(X)) stop("X is not a matrix.")
    if (!is.numeric(k)) stop("k is not numeric.")
    
    # Search predicted class(es)
    predClass <- apply(X > k, MARGIN = 1, 
                       FUN = function(row) which(row))
    
    return (predClass)
}



#' @title koplsMaxClassify
#' @description Classification function that assesses class belonging of 
#' \code{data} based on the maximum value.
#' 
#' @param X matrix. It contains the predicted response matrix Y, where 
#' columns denote classes and rows observations.
#'
#' @return
#' \item{predClass}{ matrix. It contains the predicted class(es) of \code{data}.}
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
#' 
koplsMaxClassify <- function(X) {
    # Variable format control
    if (!is.matrix(X)) stop("X is not a matrix.")
    
    # Search max position
    predClass <- apply(X = X, MARGIN = 1,
                                     FUN = function(row) which.max(row))
    
    return (predClass)
}
