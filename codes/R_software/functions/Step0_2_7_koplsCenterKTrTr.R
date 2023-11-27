#' @title koplsCenterKTrTr
#' @description Centering function for the training kernel, which is constructed
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
#' @param K matrix. Contains the training kernel matrix; K = <phi(Xtr), phi(Xtr)>. 
#'
#' @return
#' \item{K}{ matrix. The centered kernel matrix.}
#'
#' @examples
#' K <- base::matrix(stats::rnorm(n = 25), nrow = 5, ncol = 5)
#' test <- koplsCenterKTrTr_TRUE(K = K)
#' test
#' 
#' @keywords internal
#' @import stats

koplsCenterKTrTr <- function(K){
  # Variable format control
  if(!is.matrix(K)){stop("K is not a matrix.")}
  
  # Define parameters
  I <- base::diag(nrow(K))
  scaling_matrix <- base::matrix(1/nrow(K), nrow = nrow(K), ncol = ncol(K))
  
  # Center the kernel
  K <- base::crossprod(I- scaling_matrix, 
                       base::crossprod(t(K), I - scaling_matrix))
  
  # Return the centered kernel matrix
  return(K)
}
