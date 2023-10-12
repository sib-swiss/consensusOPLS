#' koplsCenterKTeTr
#' Centering function for the hybrid test/training kernel, which
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
#' @param KteTr: matrix. Contains the hybrid test/training kernel matrix, 
#' KteTr = <phi(Xte), phi(Xtr)>.
#' @param KtrTr: matrix. Contains the training kernel matrix; 
#' Ktrain = <phi(Xtr), phi(Xtr)>.
#'
#' @return
#' `KteTr`: matrix. The centered kernel matrix.
#'
#' @examples
#' KteTr <- base::matrix(1:25, nrow = 5, ncol = 5)
#' KtrTr <- base::matrix(1:25, nrow = 5, ncol = 5)
#' test <- koplsCenterKTeTr(KteTr = KteTr, KtrTr = KtrTr)
#' test

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
  
  # Calculate (1/nTrain) * I_n * I_nTrain'
  scaling_matrix <- (1/nTrain) * (I_n %*% t(I_nTrain))
  
  # Update KTeTr
  KteTr = (KteTr - scaling_matrix %*% KtrTr) %*% 
    (Itrain-(1/nTrain) * I_nTrain %*% t(I_nTrain))
  
  # Return the centered kernel matrix.
  return(KteTr)
}
