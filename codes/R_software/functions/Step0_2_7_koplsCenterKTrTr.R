#' koplsCenterKTrTr
#' Centering function for the training kernel, which is constructed
#' from the training matrix Xtr as K = <phi(Xtr), phi(Xtr)>.
#'
#' @param K: Training kernel matrix; K = <phi(Xtr), phi(Xtr)>. 
#'
#' @return
#' `K`: The centered kernel matrix.
#'
#' @examples
#' K <- matrix(1:25, nrow = 5, ncol = 5)
#' test <- koplsCenterKTrTr(K = K)
#' test
#' 
koplsCenterKTrTr <- function(K){
  # Variable format control
  if(!is.matrix(K)){stop("K is not a matrix.")}
  
  # Identity matrices
  I <- base::diag(nrow(K))
  I_n <- base::rep(x = 1, times = nrow(K))
  
  # Calculate (1/n) * I_n * I_n'
  scaling_matrix <- (1/nrow(K)) * (I_n %*% t(I_n))
  
  # Update K
  K = (I- scaling_matrix) %*% K %*% (I - scaling_matrix);
  
  # Return the centered kernel matrix
  return(K)
}
