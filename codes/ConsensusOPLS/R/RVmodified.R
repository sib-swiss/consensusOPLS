#' @title RVmodified
#' @description
#' Calculate the coefficient of determination to assess how well the OPLS 
#' regression model fits the data. It indicates how close the observed values 
#' (real data) are to the values predicted by the regression model. Its 
#' calculation is modified for the present method, as it is adjusted over the 
#' range [0, 1].
#'
#' @param X matrix. The normalized meta-kernel of all data block. 
#' @param Y matrix. The centered and/or scaled kernel of each data block.
#'
#' @return The modified R-square value.
#'
#' @examples
#' X <- matrix(runif(36), nrow=2)
#' Y <- matrix(runif(54), nrow=2)
#' result <- ConsensusOPLS:::RVmodified(X = X, Y = Y)
#' result
#' @keywords internal

RVmodified <- function(X, Y){
    # Variable format control
    if (!is.matrix(X)) stop("X is not a matrix.")
    if (!is.matrix(Y)) stop("Y is not a matrix.")
    
    AA <- tcrossprod(X)
    BB <- tcrossprod(Y)
    diag(AA) <- 0
    diag(BB) <- 0
    
    # R-square value
    RV <- sum(diag(crossprod(AA, BB))) / ((sqrt(sum(AA^2))) * (sqrt(sum(BB^2))))
    
    return(RV)
}