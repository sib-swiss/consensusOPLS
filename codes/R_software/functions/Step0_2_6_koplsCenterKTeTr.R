#' @title koplsCenterKTeTr
#' @description Centering function for the hybrid test/training kernel, which
#' is constructed from the test matrix Xte and the training matrix
#' Xtr as KteTr = <phi(Xte), phi(Xtr)>. Requires additional
#' (un-centered) training kernel to estimate mean values.
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
#' @param KteTr matrix. Contains the hybrid test/training kernel matrix, 
#' KteTr = <phi(Xte), phi(Xtr)>.
#' @param KtrTr matrix. Contains the training kernel matrix; 
#' Ktrain = <phi(Xtr), phi(Xtr)>.
#'
#' @return
#' \item{KteTr}{ matrix. The centered kernel matrix.}
#'
#' @examples
#' Xte <- base::matrix(data = stats::rnorm(n = 20), ncol=5)
#' Xtr <- base::matrix(data = stats::rnorm(n = 25), ncol=5)
#' KteTr <- koplsKernel(X1 = Xte, X2 = Xtr, Ktype='g', params=c(sigma=1.0))
#' KtrTr <- koplsKernel(X1 = Xtr, X2 = Xtr, Ktype='g', params=c(sigma=1.0))
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
  I_nTrain <- base::rep(x = 1, times = nrow(KtrTr))
  scaling_matrix <- (1/nrow(KtrTr)) * 
    base::tcrossprod(x = base::rep(x = 1, times = nrow(KteTr)), 
                     y = I_nTrain)
  
  # Center the kernel
  KteTr <- base::tcrossprod(x = KteTr - base::tcrossprod(x = scaling_matrix, 
                                                         y = t(KtrTr)),
                            y = base::diag(nrow(KtrTr)) - 1/nrow(KtrTr) * tcrossprod(I_nTrain))
  
  # Return the centered kernel matrix.
  return(KteTr)
}
