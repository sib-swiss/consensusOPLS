#' @title Variable Importance in Projection
#' @description Calculate the VIP (Variable Importance in Projection) for each
#' variable in a \code{ConsensusOPLS} model.
#' 
#' @param data A list of data blocks.
#' @param Y A vector, factor, dummy matrix or numerical matrix for the response.
#' @param model An object inheriting from class \code{ConsensusOPLS}, 
#' representing a Consensus OPLS fitted model. Default, NULL, a model will be constructed.
#' @param ... arguments to pass to \code{ConsensusOPLS}.
#'
#' @return A table with the results:
#' \code{VIP = sqrt(p*q/s)}, where
#' \code{p} is the number of variables in each block,
#' \code{q} the explained variance of Y associated to each variable, and 
#' \code{s} the total Y variance explained by the model.
#'
#' @examples
#' vip <- VIP(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")], 
#'            Y=demo_3_Omics$Y)
#' str(vip)
#' @export
#' 
VIP <- function(data, Y, model = NULL, ...) {
    # Variable format control
    if (!is.list(data)) stop("data is not a list.")
    if (!is.matrix(Y)) stop("Y is not a matrix.")
    
    # Build a model from given data if model is NULL
    if (is.null(model)) {
        rvcopls <- RVConsensusOPLS(data=data, Y=Y, ...)
        model <- rvcopls$Model
    }
    if (!is.list(model)) stop("model is not a list.")
    
    VIP <- lapply(1:length(data), function(itable) {
        # Dimensions of the data in the data
        nvariable <- ncol(data[[itable]])
        nsample   <- nrow(model$scores)
        ncomp     <- ncol(model$scores)
        
        #nclass x ncomp
        Qs <- crossprod(Y, model$scores) %*% diag(1/diag(crossprod(model$scores)), 
                                                          ncol=ncomp)
        #nsample x ncomp
        Us <- crossprod(t(Y), Qs) %*% diag(1/diag(crossprod(Qs)), 
                                           ncol=ncomp)
        #nvariable x ncomp
        Ws <- crossprod(data[[itable]], Us) %*% diag(1/diag(crossprod(Us)), 
                                                     ncol=ncomp)
        Ws <- apply(Ws, 2, function(x) x/norm(x, type='2'))
        # ncomp x ncomp x ncomp x ncomp
        s <- diag(diag(crossprod(model$scores) %*% crossprod(Qs)))
        
        VIP.itable <- t(apply(Ws, 1, function(x) {
            q <- crossprod(s, x^2)
            return (sqrt(nvariable * q / sum(s)))
        }))
        colnames(VIP.itable) <- colnames(model$scores)
        
        return (VIP.itable)
    })
    names(VIP) <- names(data)
    
    return (VIP)
}
