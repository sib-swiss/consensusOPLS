#' koplsCenterKTeTr
#' Centering function for the hybrid test/training kernel, which
#' is constructed from the test matrix Xte and the training matrix
#' Xtr as KteTr = <phi(Xte), phi(Xtr)>. Requires additional
#' (un-centered) training kernel to estimate mean values.
#'
#' @param KteTr: hybrid test/training kernel matrix, KteTr = <phi(Xte), phi(Xtr)>.
#' @param KtrTr: Training kernel matrix; Ktrain = <phi(Xtr), phi(Xtr)>.
#'
#' @return
#' `KteTr`: The centered kernel matrix.
#'
#' @examples
#' KteTr <- matrix(1:25, nrow = 5, ncol = 5)
#' KtrTr <- matrix(1:25, nrow = 5, ncol = 5)
#' test <- koplsCenterKTeTr(KteTr = KteTr, KtrTr = KtrTr)
#' test

koplsCenterKTeTr <- function(KteTr, KtrTr){
  # Variable format control
  if(!is.matrix(KteTr)){stop("KteTr is not a matrix.")}
  if(!is.matrix(KtrTr)){stop("KtrTr is not a matrix.")}
  
  # Define parameters
  Itrain <-  base::diag(ncol(KtrTr))
  I_nTrain <- base::rep(x = 1, times = ncol(KtrTr))
  nTrain <- ncol(KtrTr)
  
  I <- base::diag(ncol(KteTr))
  I_n <- base::rep(x = 1, times = ncol(KteTr))
  n <- ncol(KteTr)
  
  # Calculate (1/nTrain) * I_n * I_nTrain'
  scaling_matrix <- (1/nTrain) * (I_n %*% t(I_nTrain))
  
  # Update KTeTr
  KteTr = (KteTr - scaling_matrix %*% KtrTr) %*% 
    (Itrain-(1/nTrain) * I_nTrain %*% t(I_nTrain))
  
  # Return the centered kernel matrix.
  return(KteTr)
}
