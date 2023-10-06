#' RVmodified
#' Calculates the coefficient of determination to assess how well the OPLS 
#' regression model fits the data. It indicates how close the observed values 
#' (real data) are to the values predicted by the regression model. Its 
#' calculation is modified for the present method, as it is adjusted over the 
#' range [0, 1].
#'
#' @param X: the normalized meta-kernel of all data block. 
#' @param Y: the centered and/or scaled kernel of each data block.
#'
#' @return
#' `RV`: the modified R-square value
#'
#' @examples
#' TO DO

RVmodified <- function(X, Y){
  # Variable format control
  if(!is.matrix(X)){stop("X is not a matrix.")}
  if(!is.matrix(Y)){stop("Y is not a matrix.")}
  
  # Dimension reduction
  AA <- X%*%t(X)
  BB <- Y%*%t(Y)
  
  # Similarity matrix
  AA0 <- AA - diag(diag(AA), nrow(AA), ncol(AA))
  BB0 <- BB - diag(diag(BB), nrow(BB), ncol(BB))
  
  # Return the R-square value
  RV <- sum(AA0 %*% BB0) / (sqrt(sum(AA0**2)) * sqrt(sum(BB0**2)))
  return(RV)
}
