#' @title koplsCenterKTeTe
#' @description Centering function for the test kernel, which is constructed
#' from the test matrix Xte as KteTe = <phi(Xte), phi(Xte)>.
#' Requires additional (un-centered) kernels KteTr and KteTr to estimate mean 
#' values.
#' 
#' # ------------------------------------------------------------------------ #
#' This file is part of the K-OPLS package, developed by Max Bylesjo, 
#' University of Umea, Judy Fonville and Mattias Rantalainen, Imperial College.
#' 
#' Copyright (c) 2007-2010 Max Bylesjo, Judy Fonville and Mattias Rantalainen 
#' 
#' This code has been extended and adapted under the terms of the GNU General 
#' Public License version 2 as published by the Free Software Foundation.
#' # ------------------------------------------------------------------------ #
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
#' KteTr <- matrix(1:25, nrow = 5, ncol = 5)
#' KteTe <- matrix(1:25, nrow = 5, ncol = 5)
#' KtrTr <- matrix(1:25, nrow = 5, ncol = 5)
#' test <- koplsCenterKTeTe(KteTe = KteTe, KteTr = KteTr, KtrTr = KtrTr)
#' test
#' 
#' @keywords internal

koplsCenterKTeTe <- function(KteTe, KteTr, KtrTr){
  # Variable format control
  if (!is.matrix(KteTe) || !is.matrix(KteTr) || !is.matrix(KtrTr)) {
    stop("One or more inputs are not matrices.")
  }
  
  # Define parameters
  scaling_matrix <- (1/ncol(KteTr)) * rep(x = 1, 
                                          times = nrow(KteTr)) %*% 
    t(rep(x = 1, times = ncol(KteTr)))
  
  # Center the kernel
  KteTe <- KteTe - scaling_matrix %*% t(KteTr) - KteTr %*% t(scaling_matrix) + 
    scaling_matrix %*% KtrTr %*% t(scaling_matrix)
  
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
#' KteTr <- matrix(1:25, nrow = 5, ncol = 5)
#' KtrTr <- matrix(1:25, nrow = 5, ncol = 5)
#' test <- koplsCenterKTeTr(KteTr = KteTr, KtrTr = KtrTr)
#' test
#' 
#' @keywords internal

koplsCenterKTeTr <- function(KteTr, KtrTr){
  # Variable format control
  if (!is.matrix(KteTr) || !is.matrix(KtrTr)) {
    stop("One or more inputs are not matrices.")
  }
  
  # Define parameters
  I_nTrain <- rep(x = 1, times = nrow(KtrTr))
  scaling_matrix <- (1/nrow(KtrTr)) * rep(x = 1, 
                                          times = nrow(KteTr)) %*% t(I_nTrain)
  
  # Center the kernel
  KteTr <- (KteTr - scaling_matrix %*% KtrTr) %*% 
    (diag(nrow(KtrTr)) - 1/nrow(KtrTr) * I_nTrain %*% t(I_nTrain))
  
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
#' test <- koplsCenterKTrTr(K = K)
#' test
#' 
#' @keywords internal

koplsCenterKTrTr <- function(K){
  # Variable format control
  if(!is.matrix(K)){stop("K is not a matrix.")}
  
  # Define parameters
  I <- diag(nrow(K))
  scaling_matrix <- (1/nrow(K)) * (rep(x = 1, times = nrow(K)) %*% 
                                     t(rep(x = 1, times = nrow(K))))
  
  # Center the kernel
  K <- (I- scaling_matrix) %*% K %*% (I - scaling_matrix)
  
  # Return the centered kernel matrix
  return(K)
}
