#' RVmodified
#' Calculates the coefficient of determination to assess how well the OPLS 
#' regression model fits the data. It indicates how close the observed values 
#' (real data) are to the values predicted by the regression model. Its 
#' calculation is modified for the present method, as it is adjusted over the 
#' range [0, 1].
#'
#' @param X: matrix. The normalized meta-kernel of all data block. 
#' @param Y: matrix. The centered and/or scaled kernel of each data block.
#'
#' @return
#' `RV`: numeric. The modified R-square value.
#'
#' @examples
#' X <- base::matrix(c(1, 2, 3, 4, 5, 6), nrow = 2)
#' Y <- base::matrix(c(2, 4, 6, 8, 10, 12), nrow = 2)
#' result <- RVmodified(X = X, Y = Y)
#' cat("R-Square Value:", result, "\n")

RVmodified <- function(X, Y){
  # Variable format control
  if(!is.matrix(X)){stop("X is not a matrix.")}
  if(!is.matrix(Y)){stop("Y is not a matrix.")}
  
  # Dimension reduction
  AA <- X %*% t(X)
  BB <- Y %*% t(Y)
  
  # Similarity matrix
  AA0 <- AA - base::diag( base::diag(AA), nrow(AA), ncol(AA))
  BB0 <- BB - base::diag( base::diag(BB), nrow(BB), ncol(BB))
  
  # Return the R-square value
  RV <- base::sum( base::diag(AA0 %*% BB0)) / 
    (base::sqrt( base::sum(AA0**2))) / 
    (base::sqrt( base::sum(BB0**2)))
  return(RV)
}
