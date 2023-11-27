#' @title koplsCenterKTeTe
#' @description Centering function for the test kernel, which is constructed
#' from the test matrix Xte as KteTe = <phi(Xte), phi(Xte)>.
#' Requires additional (un-centered) kernels KteTr and KteTr to estimate mean 
#' values.
#'
#' @param KteTe matrix. Contains the test kernel matrix, 
#' KteTe = <phi(Xte), phi(Xte)>.  
#' @param KteTr matrix. Contains the hybrid test/training kernel matrix, 
#' KteTr = <phi(Xte), phi(Xtr)>.
#' @param KtrTr matrix. Contains the training kernel matrix, 
#' KtrTr = <phi(Xtr), phi(Xtr)>.
#'
#' @return
#' \item{KteTe}{ matrix. The centered test kernel matrix.}
#'
#' @examples
#' KteTr <- matrix(rnorm(n = 25), nrow = 5, ncol = 5)
#' KteTe <- matrix(rnorm(n = 25), nrow = 5, ncol = 5)
#' KtrTr <- matrix(rnorm(n = 25), nrow = 5, ncol = 5)
#' test <- ConsensusOPLS:::koplsCenterKTeTe(KteTe = KteTe, 
#'                                          KteTr = KteTr, 
#'                                          KtrTr = KtrTr)
#' test
#' 
#' @keywords internal
#' @import stats

koplsCenterKTeTe <- function(KteTe, KteTr, KtrTr){
  # Variable format control
  if (!is.matrix(KteTe) || !is.matrix(KteTr) || !is.matrix(KtrTr)) {
    stop("One or more inputs are not matrices.")
  }
  
  # Define parameters
  scaling_matrix <- matrix(1/nrow(KteTr), nrow = nrow(KteTr), ncol = ncol(KteTr))
  
  # Center the kernel
  KteTe <- KteTe - tcrossprod(scaling_matrix, KteTr) - 
    tcrossprod(KteTr, scaling_matrix) + 
    tcrossprod(scaling_matrix, tcrossprod(scaling_matrix, KtrTr))
  
  # Return the centered test kernel matrix
  return(KteTe)
}



#' @title koplsCenterKTeTr
#' @description Centering function for the hybrid test/training kernel, which
#' is constructed from the test matrix Xte and the training matrix
#' Xtr as KteTr = <phi(Xte), phi(Xtr)>. Requires additional
#' (un-centered) training kernel to estimate mean values.
#'
#' @param KteTr matrix. Contains the hybrid test/training kernel matrix, 
#' KteTr = <phi(Xte), phi(Xtr)>.
#' @param KtrTr matrix. Contains the training kernel matrix; 
#' Ktrain = <phi(Xtr), phi(Xtr)>.
#' 
#' @return
#' \item{KteTr}{ matrix. The centered kernel matrix.}
#' 
#' @examples
#' KteTr <- matrix(rnorm(n = 25), nrow = 5, ncol = 5)
#' KtrTr <- matrix(rnorm(n = 25), nrow = 5, ncol = 5)
#' test <- ConsensusOPLS:::koplsCenterKTeTr(KteTr = KteTr, KtrTr = KtrTr)
#' test
#' 
#' @keywords internal
#' @import stats

koplsCenterKTeTr <- function(KteTr, KtrTr){
  # Variable format control
  if (!is.matrix(KteTr) || !is.matrix(KtrTr)) {
    stop("One or more inputs are not matrices.")
  }
  
  # Define parameters
  I_nTrain <- rep(x = 1, times = nrow(KtrTr))
  scaling_matrix <- (1/nrow(KtrTr)) * tcrossprod(rep(x = 1, 
                                                     times = nrow(KteTr)), 
                                                 I_nTrain)
  
  # Center the kernel
  KteTr <- tcrossprod(KteTr - tcrossprod(scaling_matrix, t(KtrTr)),
                      diag(nrow(KtrTr)) - 1/nrow(KtrTr) * tcrossprod(I_nTrain))
  
  # Return the centered kernel matrix.
  return(KteTr)
}



#' @title koplsCenterKTrTr
#' @description Centering function for the training kernel, which is constructed
#' from the training matrix Xtr as K = <phi(Xtr), phi(Xtr)>.
#'
#' @param K matrix. Contains the training kernel matrix; K = <phi(Xtr), phi(Xtr)>. 
#'
#' @return
#' \item{K}{ matrix. The centered kernel matrix.}
#'
#' @examples
#' K <- matrix(1:25, nrow = 5, ncol = 5)
#' test <- ConsensusOPLS:::koplsCenterKTrTr(K = K)
#' test
#' 
#' @keywords internal
#' @import stats

koplsCenterKTrTr <- function(K){
  # Variable format control
  if(!is.matrix(K)){stop("K is not a matrix.")}
  
  # Define parameters
  I <- diag(nrow(K))
  scaling_matrix <- matrix(1/nrow(K), nrow = nrow(K), ncol = ncol(K))
  
  # Center the kernel
  K <- crossprod(I- scaling_matrix, crossprod(t(K), I - scaling_matrix))
  
  # Return the centered kernel matrix
  return(K)
}
