#' koplsCenterKTrTr
#' Centering function for the training kernel, which is constructed
#' from the training matrix Xtr as K = <phi(Xtr), phi(Xtr)>.
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
#' @param K: matrix. Contains the training kernel matrix; K = <phi(Xtr), phi(Xtr)>. 
#'
#' @return
#' `K`: matrix. The centered kernel matrix.
#'
#' @examples
#' K <- base::matrix(1:25, nrow = 5, ncol = 5)
#' test <- koplsCenterKTrTr(K = K)
#' test

koplsCenterKTrTr <- function(K){
  # Variable format control
  if(!is.matrix(K)){stop("K is not a matrix.")}
  
  # Identity matrices
  I <- base::diag(nrow(K))
  I_n <- base::rep(x = 1, times = nrow(K))
  
  # Calculate (1/n) * I_n * I_n'
  scaling_matrix <- (1/nrow(K)) * (I_n %*% t(I_n))
  
  # Update K
  K = (I- scaling_matrix) %*% K %*% (I - scaling_matrix)
  
  # Return the centered kernel matrix
  return(K)
}
