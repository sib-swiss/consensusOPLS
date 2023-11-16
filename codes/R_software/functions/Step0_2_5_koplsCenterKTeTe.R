#' koplsCenterKTeTe
#' Centering function for the test kernel, which is constructed
#' from the test matrix Xte as KteTe = <phi(Xte), phi(Xte)>.
#' Requires additional (un-centered) kernels KteTr and KteTr to
#' estimate mean values.
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
#' @param KteTe: matrix. Contains the test kernel matrix, 
#' KteTe = <phi(Xte), phi(Xte)>.  
#' @param KteTr: matrix. Contains the hybrid test/training kernel matrix, 
#' KteTr = <phi(Xte), phi(Xtr)>.
#' @param KtrTr: matrix. Contains the training kernel matrix, 
#' KtrTr = <phi(Xtr), phi(Xtr)>.
#'
#' @return
#' `KteTe`: matrix. The centered test kernel matrix.
#'
#' @examples
#' KteTr <- base::matrix(1:25, nrow = 5, ncol = 5)
#' KteTe <- matrix(1:25, nrow = 5, ncol = 5)
#' KtrTr <- matrix(1:25, nrow = 5, ncol = 5)
#' test <- koplsCenterKTeTe(KteTe = KteTe, KteTr = KteTr, KtrTr = KtrTr)

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
  D_te <- (1/nTrain) * I_n %*% t(I_nTrain)
  KteTe <- KteTe - D_te%*%t(KteTr) - KteTr%*%t(D_te) + D_te%*%KtrTr%*%t(D_te)
  
  # Return the centered test kernel matrix
  return(KteTe)
}
