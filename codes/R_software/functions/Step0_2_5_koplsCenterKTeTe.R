#' koplsCenterKTeTe
#' Centering function for the test kernel, which is constructed
#' from the test matrix Xte as KteTe = <phi(Xte), phi(Xte)>.
#' Requires additional (un-centered) kernels KteTr and KteTr to
#' estimate mean values.
#'
#' @param KteTe: Test kernel matrix, KteTe = <phi(Xte), phi(Xte)>.  
#' @param KteTr: Test/training kernel matrix, KteTr = <phi(Xte), phi(Xtr)>.
#' @param KtrTr: Training kernel matrix, KtrTr = <phi(Xtr), phi(Xtr)>.
#'
#' @return
#' `KteTe` the centered test kernel matrix
#'
#' @examples
#' # TO DO

koplsCenterKTeTe <- function(KteTe, KteTr, KtrTr){
  # Variable format control
  if(!is.matrix(KteTe)){stop("KteTe is not a matrix.")}
  if(!is.matrix(KteTr)){stop("KteTr is not a matrix.")}
  if(!is.matrix(KtrTr)){stop("KtrTr is not a matrix.")}
  
  # Define parameters
  Itrain <- base::diag(ncol(KteTr))
  I_nTrain <- base::rep(x = 1, times = ncol(KteTr))
  nTrain <- ncol(KteTr)
  
  I <- base::diag(nrow(KteTr))
  I_n <- base::rep(x = 1, times = nrow(KteTr))
  n <- nrow(KteTr)
  
  # Center the kernel
  D_te = (1/nTrain) %.*% I_n %*% t(I_nTrain)
  KteTe = KteTe - D_te%*%t(KteTr) - KteTr%*%t(D_te) + D_te%*%KtrTr%*%t(D_te)
  
  # Return the centered test kernel matrix
  return(KteTe)
}