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
#' # TO DO

koplsCenterKTeTr <- function(KteTr, KtrTr){
  # Variable format control
  if(!is.matrix(KteTr)){stop("KteTr is not a matrix.")}
  if(!is.matrix(KtrTr)){stop("KtrTr is not a matrix.")}
  
  # Define parameters
  Itrain <-  base::diag(nrow(KtrTr))
  I_nTrain <- base::rep(x = 1, times = nrow(KtrTr))
  nTrain <- nrow(KtrTr)
  
  I <- base::diag(nrow(KteTr))
  I_n <- base::rep(x = 1, times = nrow(KteTr))
  n <- nrow(KteTr)
  
  KteTr = (KteTr-(1/nTrain)%*% I_n %*% t(I_nTrain) %*% KtrTr) %*% 
    (Itrain-(1/nTrain) %.*% I_nTrain %*% t(I_nTrain))
  
  # Return the centered kernel matrix.
  return(KteTr)
}