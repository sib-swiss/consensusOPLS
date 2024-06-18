#' @title koplsKernel
#' @description 
#' This function constructs a kernel matrix \code{K = <phi(X1), phi(X2)>}.
#' The kernel function \code{phi} determines how the data is transformed and
#' is passed as the separate parameter \code{type} to the function.
#' Currently \code{type} can be either \code{g} for Gaussian or \code{p} for
#' polynomial.
#'
#' @param X1 matrix. The first X matrix (non-centered). This is the left side in 
#' the expression K = <phi(X1), phi(X2)>.
#' @param X2 matrix. The second X matrix (non-centered). This is the right side
#' in the expression K = <phi(X1), phi(X2)>. If X2 = [] (empty set), then only
#' X1 will be used for the calculations. This way, only (n^2 - n)/2 instead of
#' n^2 calculations have to be performed, which is typically much faster. Only
#' applicable for pure training or testing kernels.
#' @param type character. Indicates the type of kernel used. Supported entries 
#' are \code{g} for Gaussian kernel, and \code{p} for Polynomial kernel. 
#' @param params vector. It contains parameter for the kernel function. 
#' Currently, all supported kernel functions use a scalar value so the vector 
#' property of the parameters is for future compatibility.
#' 
#' @returns The kernel matrix transformed by the kernel function specified by 
#' \code{type}.
#'
#' @examples
#' X1 <- base::matrix(stats::rnorm(n = 20), nrow = 5)
#' X2 <- base::matrix(stats::rnorm(n = 24), nrow = 6)
#' 
#' # Polynomial example
#' params_polynomial <- c(order=2)  # Polynomial kernel order
#' kernel_polynomial <- ConsensusOPLS:::koplsKernel(X1 = X1, X2 = X2, 
#'                                                  type = 'p', params=params_polynomial)
#' kernel_polynomial
#' 
#' @keywords internal
#' @noRd
#' 
koplsKernel <- function(X1, X2 = NULL, type = 'p', params = c(order=1.0)) {
    # Variable format control
    if (!is.matrix(X1)) 
        stop("X1 is not a matrix.")
    if (!is.null(X2) && !is.matrix(X2)) 
        stop("X2 is not a matrix.")
    if (!is.vector(params) || !is.numeric(params)) 
        stop("params must be a numeric vector.")
    type <- match.arg(type, choices=c('p', 'g'))
    
    if (type == "g") { # Define Gaussian Kernel
        #TODO: check why the kernel matrix becomes identity matrix when there are more variables than samples
        sigma <- params['sigma'] #small value = overfit, larger = more general
        
        if (is.null(X2) || nrow(X2) == 0) {
            K <- exp(-as.matrix(dist(X1, method = "euclidean"))^2/(2*(sigma^2)))
        } else{
            K <- exp(-as.matrix(dist(rbind(X1, X2), 
                                     method = "euclidean"))[1:nrow(X1), 
                                                            nrow(X1)+1:nrow(X2)]^2/(2*(sigma^2)))
        }
    } else if (type == "p") { # Define Polynomial Kernel
        porder <- params['order']
        
        if (is.null(X2) || nrow(X2) == 0) {
            K <- (tcrossprod(X1) + 1)^porder
        } else {
            K <- (tcrossprod(x = X1, y = X2) + 1)^porder
        }
    }
    
    # Return the kernel matrix
    return (K)
}
