#' koplsKernel
#' Constructs a kernel matrix K = <phi(X1), phi(X2)>.
#' The kernel function k() determines how the data is transformed and
#' is passed as the separate parameter 'Ktype' to the function.
#' Currently 'Ktype' can be either 'g' (Gaussian) or 'p' (polynomial).
#'
#' @param X1: the first X matrix (non-centered). This is the left side in the 
#' expression K = <phi(X1), phi(X2)>.
#' @param X2: the second X matrix (non-centered). This is the right side in the 
#' expression K = <phi(X1), phi(X2)>. If X2 = [] (empty set), then only X1 will 
#' be used for the calculations. This way, only (n^2 - n)/2 instead of n^2 
#' calculations have to be performed, which is typically much faster. Only 
#' applicable for pure training or testing kernels.
#' @param Ktype: the type of kernel used. Supported entries are `g` for Gaussian 
#' kernel, and `p` for Polynomial kernel. 
#' @param params: A vector with parameter for the kernel function. (Currently, 
#' all supported kernel functions use a scalar value so the vector property of 
#' the parameters is for future compability).
#' 
#' @return
#' `K`: The kernel matrix transformed by the kernel function specified by `Ktype`.
#'
#' @examples
#' TO DO

koplsKernel <- function(X1, X2 = NULL, Ktype = 'g', params){
  # Variable format control
  if(!is.matrix(X1)){stop("X1 is not a matrix.")}
  if(!is.null(X2)){if(!is.matrix(X2)){stop("X2 is not a matrix.")}}
  if(!is.character(Ktype)){
    stop("Ktype is not a character.")
  } else{
    if(!(Ktype %in% c('g', 'p'))){
      stop("Ktype must be `g` or `p`.")
    }
  }
  if(!is.vector(params)){stop("params must be a vector.")}
  
  # Test if there is only one matrix for initialize the kernel matrix
  if(nrow(X2) == 0 | is.null(X2)){
    K <- matrix(0, nrow = nrow(X1), ncol = nrow(X1))
  } else {
    K <- matrix(0, nrow = nrow(X1), ncol = nrow(X2))
  }
  
  # Define Gaussian Kernel
  if(Ktype == "g"){
    sigma <- params[1] #small value = overfit, larger = more general
    
    if(nrow(X2) == 0 | is.null(X2)){
      # Calculate upper triangle and duplicated due to symmetry
      for(i in 1:nrow(X1)){
        for(j in 1:nrow(X1)){
          diff_norm <- sum((X1[i,] - X1[j,])**2)
          K[i,j] <- exp(-diff_norm/ (2*(sigma**2)) )
          K[j, i] <- K[i, j]
        }
      }
    } else{
      # Loop over the entire kernel matrix
      for (i in 1:nrow(X1)) {
        for (j in 1:nrow(X2)) {
          diff_norm <- sum((X1[i,] - X2[j,])**2)
          K[i, j] <- exp(-diff_norm / (2 * (sigma**2)))
        }
      }
    }
  }
  
  # Define Polnomial Kernel
  if(Ktype == "p"){
    porder <- params[1]
    
    if (nrow(X2) == 0 | is.null(X2)) {
      # Calculate upper triangle and duplicate due to symmetry
      for (i in 1:nrow(X1)) {
        for (j in i:nrow(X1)) {
          K[i, j] <- (sum(X1[i,] * X1[j,]) + 1)**porder
          K[j, i] <- K[i, j]
        }
      }
    } else {
      # Loop over the entire kernel matrix
      for (i in 1:nrow(X1)) {
        for (j in 1:nrow(X2)) {
          K[i, j] <- (sum(X1[i,] * X2[j,]) + 1)**porder
        }
      }
    }
  }
  
  # Return the kernel matrix transformed
  return(K)
}
