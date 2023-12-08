#' @title koplsKernel
#' @description 
#' This function constructs a kernel matrix \code{K = <phi(X1), phi(X2)>}.
#' The kernel function \code{phi} determines how the data is transformed and
#' is passed as the separate parameter 'Ktype' to the function.
#' Currently 'Ktype' can be either 'g' (Gaussian) or 'p' (polynomial).
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
#' @param X1 matrix. The first X matrix (non-centered). This is the left side in the 
#' expression K = <phi(X1), phi(X2)>.
#' @param X2 matrix. The second X matrix (non-centered). This is the right side in the 
#' expression K = <phi(X1), phi(X2)>. If X2 = [] (empty set), then only X1 will 
#' be used for the calculations. This way, only (n^2 - n)/2 instead of n^2 
#' calculations have to be performed, which is typically much faster. Only 
#' applicable for pure training or testing kernels.
#' @param Ktype character. Indicates the type of kernel used. Supported entries 
#' are \code{g} for Gaussian kernel, and \code{p} for Polynomial kernel. 
#' @param params vector. It contains parameter for the kernel function. (Currently, 
#' all supported kernel functions use a scalar value so the vector property of 
#' the parameters is for future compatibility).
#' 
#' @return The kernel matrix transformed by the kernel function specified by 
#' \code{Ktype}.
#'
#' @examples
#' X1 <- base::matrix(stats::rnorm(n = 20), nrow = 5)
#' X2 <- base::matrix(stats::rnorm(n = 24), nrow = 6)
#' # Gaussian example
#' params_gaussian <- c(sigma=1.0)  # Sigma values
#' kernel_gaussian <- koplsKernel(X1, X2, Ktype='g', params=params_gaussian)
#' kernel_gaussian
#' 
#' # Polynomial example
#' params_polynomial <- c(order=2)  # Polynomial kernel order
#' kernel_polynomial <- koplsKernel(X1, X2, Ktype='p', params=params_polynomial)
#' kernel_polynomial
#' 
#' @keywords internal

koplsKernel <- function(X1, X2 = NULL, Ktype = 'g', params = c(sigma=1.0)){
  # Variable format control
  if (!is.matrix(X1)) 
    stop("X1 is not a matrix.")
  if (!is.null(X2) && !is.matrix(X2)) 
    stop("X2 is not a matrix.")
  if (!is.vector(params) || !is.numeric(params)) 
    stop("params must be a numeric vector.")
  Ktype <- match.arg(Ktype, choices=c('g', 'p'))
  
  # Test if there is only one matrix for initialize the kernel matrix
  if( isTRUE(nrow(X2) == 0) | is.null(X2)){
    K <- base::matrix(data = 0, nrow = nrow(X1), ncol = nrow(X1))
  } else {
    K <- base::matrix(data = 0, nrow = nrow(X1), ncol = nrow(X2))
  }
  
  if (Ktype == "g") { # Define Gaussian Kernel
    sigma <- params['sigma'] #small value = overfit, larger = more general
    
    if (is.null(X2) || nrow(X2) == 0) {
      K <- exp(-as.matrix(dist(X1, method = "euclidean"))^2/(2*(sigma^2)))
    } else{
      K <- exp(-as.matrix(dist(rbind(X1, X2), 
                               method = "euclidean"))[1:nrow(X1), 
                                                      nrow(X1)+1:nrow(X2)]^2/(2*(sigma^2)))
    }
  } else if (Ktype == "p") { # Define Polnomial Kernel
    porder <- params['order']
    
    if (is.null(X2) || nrow(X2) == 0) {
      K <- (base::tcrossprod(X1) + 1)^porder
    } else {
      K <- (base::tcrossprod(X1, X2) + 1)^porder
    }
  }
  
  # Return the kernel matrix transformed
  return(K)
}
