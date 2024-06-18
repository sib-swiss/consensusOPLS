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
#' @returns The modified R-square value.
#'
#' @examples
#' X <- matrix(stats::rnorm(n = 3600), nrow=200)
#' Y <- matrix(stats::rnorm(n = 5400), nrow=200)
#' result <- ConsensusOPLS:::RVmodified(X = X, Y = Y)
#' result
#' 
#' @keywords internal
#' @noRd
#' 
RVmodified <- function(X, Y){
    # Variable format control
    if (!is.matrix(X)) stop("X is not a matrix.")
    if (!is.matrix(Y)) stop("Y is not a matrix.")
    
    AA <- tcrossprod(X)
    BB <- tcrossprod(Y)
    diag(AA) <- 0
    diag(BB) <- 0
    
    # R-square value
    RV <- sum(diag(crossprod(x = AA, y = BB))) / 
        ((sqrt(sum(AA^2))) * (sqrt(sum(BB^2))))
    
    return (RV)
}
