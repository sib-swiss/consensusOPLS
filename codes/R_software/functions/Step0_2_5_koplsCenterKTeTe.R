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
#' Xte <- base::matrix(data = stats::rnorm(n = 20), ncol=5)
#' Xtr <- base::matrix(data = stats::rnorm(n = 25), ncol=5)
#' KteTe <- koplsKernel(X1 = Xte, X2 = Xte, Ktype='g', params=c(sigma=1.0))
#' KteTr <- koplsKernel(X1 = Xte, X2 = Xtr, Ktype='g', params=c(sigma=1.0))
#' KtrTr <- koplsKernel(X1 = Xtr, X2 = Xtr, Ktype='g', params=c(sigma=1.0))
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
  scaling_matrix <- base::matrix(data = 1/nrow(KteTr), nrow = nrow(KteTr), 
                                 ncol = ncol(KteTr))
  
  # Center the kernel
  KteTe <- KteTe - base::tcrossprod(x = scaling_matrix, y = KteTr) - 
    base::tcrossprod(x = KteTr, y = scaling_matrix) + 
    base::tcrossprod(x = scaling_matrix, y = base::tcrossprod(x = scaling_matrix, 
                                                              y = KtrTr))
  
  # Return the centered test kernel matrix
  return(KteTe)
}
