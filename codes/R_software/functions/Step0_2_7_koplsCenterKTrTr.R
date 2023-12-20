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
#' Xtr <- base::matrix(data = stats::rnorm(n = 25), ncol=5)
#' K <- koplsKernel(X1 = Xtr, X2 = Xtr, Ktype='g', params=c(sigma=1.0))
#' test <- koplsCenterKTrTr_TRUE(K = K)
#' test
#' 
#' @keywords internal

koplsCenterKTrTr <- function(K){
  # Variable format control
  if(!is.matrix(K)){stop("K is not a matrix.")}
  
  # Define parameters
  scaling_matrix <- base::diag(nrow(K)) - 1/nrow(K)
  
  # Center the kernel
  K <- base::crossprod(x = scaling_matrix, 
                       y = base::crossprod(x = t(K), y = scaling_matrix))
  
  # Return the centered kernel matrix
  return(K)
}
