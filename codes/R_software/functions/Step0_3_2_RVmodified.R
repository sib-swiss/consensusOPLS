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
#' X <- base::matrix(stats::runif(n = 36), nrow = 2)
#' Y <- base::matrix(stats::runif(n = 54), nrow = 2)
#' result <- RVmodified(X = X, Y = Y)
#' result
#' 
#' @keywords internal

RVmodified <- function(X, Y){
  # Variable format control
  if(!is.matrix(X)){stop("X is not a matrix.")}
  if(!is.matrix(Y)){stop("Y is not a matrix.")}
  
  # Dimension reduction
  AA <- tcrossprod(X)
  BB <- tcrossprod(Y)
  
  # Similarity matrix
  diag(AA) <- 0
  diag(BB) <- 0
  
  # Return the R-square value
  # R-square value
  RV <- sum(diag(crossprod(x = AA, y = BB))) / 
    ((sqrt(sum(AA^2))) * (sqrt(sum(BB^2))))

  # Return the result
  return(RV)
}
