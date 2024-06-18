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
#' @returns The centered test kernel matrix.
#' @import stats
#' @examples
#' Xte <- matrix(data = stats::rnorm(n = 20), ncol=5)
#' Xtr <- matrix(data = stats::rnorm(n = 25), ncol=5)
#' KteTe <- ConsensusOPLS:::koplsKernel(X1 = Xte, X2 = Xte, 
#'                                      type='p', params=c(order=1.0))
#' KteTr <- ConsensusOPLS:::koplsKernel(X1 = Xte, X2 = Xtr, 
#'                                      type='p', params=c(order=1.0))
#' KtrTr <- ConsensusOPLS:::koplsKernel(X1 = Xtr, X2 = Xtr, 
#'                                      type='p', params=c(order=1.0))
#' ConsensusOPLS:::koplsCenterKTeTe(KteTe = KteTe, 
#'                                  KteTr = KteTr, 
#'                                  KtrTr = KtrTr)
#' @keywords internal
#' @noRd
#' 
koplsCenterKTeTe <- function(KteTe, KteTr, KtrTr) {
    # Variable format control
    if (!is.matrix(KteTe) || !is.matrix(KteTr) || !is.matrix(KtrTr)) {
        stop("One or more inputs are not matrices.")
    }
    
    # Define parameters
    scaling_matrix <- matrix(data = 1/nrow(KteTr), nrow = nrow(KteTr), 
                             ncol = ncol(KteTr))
    
    # Center the kernel
    KteTe <- KteTe - 
        tcrossprod(x = scaling_matrix, y = KteTr) - 
        tcrossprod(x = KteTr, y = scaling_matrix) + 
        tcrossprod(x = scaling_matrix, y = tcrossprod(x = scaling_matrix, 
                                                      y = KtrTr))
    
    # Return the centered test kernel matrix
    return (KteTe)
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
#' @returns The centered kernel matrix.
#' 
#' @examples
#' Xte <- matrix(data = stats::rnorm(n = 20), ncol=5)
#' Xtr <- matrix(data = stats::rnorm(n = 25), ncol=5)
#' KteTr <- ConsensusOPLS:::koplsKernel(X1 = Xte, X2 = Xtr, 
#'                                      type='p', params=c(order=1.0))
#' KtrTr <- ConsensusOPLS:::koplsKernel(X1 = Xtr, X2 = Xtr, 
#'                                      type='p', params=c(order=1.0))
#' ConsensusOPLS:::koplsCenterKTeTr(KteTr = KteTr, KtrTr = KtrTr)
#' 
#' @keywords internal
#' @noRd
#' 
koplsCenterKTeTr <- function(KteTr, KtrTr) {
    # Variable format control
    if (!is.matrix(KteTr) || !is.matrix(KtrTr) || !is.numeric(KteTr) || !is.numeric(KtrTr)) {
        stop("One or more inputs are not numeric matrices.")
    }
    
    # Define parameters
    I_nTrain <- rep(x = 1, times = nrow(KtrTr))
    scaling_matrix <- (1/nrow(KtrTr)) * tcrossprod(x = rep(x = 1, 
                                                           times = nrow(KteTr)), 
                                                   y = I_nTrain)
    
    # Center the kernel
    KteTr <- tcrossprod(x = KteTr - tcrossprod(x = scaling_matrix, 
                                               y = t(KtrTr)),
                        y = diag(nrow(KtrTr)) - 1/nrow(KtrTr) * tcrossprod(I_nTrain))
    
    # Return the centered kernel matrix.
    return (KteTr)
}



#' @title koplsCenterKTrTr
#' @description Centering function for the training kernel, which is constructed
#' from the training matrix Xtr as K = <phi(Xtr), phi(Xtr)>.
#'
#' @param K matrix. Contains the training kernel matrix; K = <phi(Xtr), phi(Xtr)>. 
#'
#' @returns The centered kernel matrix.
#' @import stats
#' @examples
#' Xtr <- matrix(data = stats::rnorm(n = 25), ncol = 5)
#' K <- ConsensusOPLS:::koplsKernel(X1 = Xtr, X2 = Xtr, 
#'                                  type='p', params=c(order=1.0))
#' ConsensusOPLS:::koplsCenterKTrTr(K = K)
#' 
#' @keywords internal
#' @noRd
#' 
koplsCenterKTrTr <- function(K) {
    return (scale(t(scale(K, scale=F)), scale=F))
}
