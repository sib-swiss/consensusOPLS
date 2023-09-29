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
#' # TO DO
#' 
koplsCenterKTrTr <- function(K){
  I <- base::diag(nrow(K))
  I_n <- base::rep(x = 1, times = nrow(K))
  K = (I- (1/nrow(K)) %.*% I_n %*% t(I_n)) %*% K %*% 
    (I-(1/nrow(K)) %.*% I_n %*% t(I_n));
  
  # Return the centered kernel matrix
  return(K)
}