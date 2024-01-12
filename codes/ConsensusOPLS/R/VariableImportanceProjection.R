#' @title MBVIP
#' @description Calculate the VIP (Variable Importance in Projection) for each variable of a 
#' Consensus-OPLS model.
#' 
#' @param data The collection list containing each block of data.
#' @param Y matrix. The response matrix to predict.
#' @param model A consensusOPLS model. Default, NULL, the model will be built.
#' @param mc.cores Number of cores for parallel computing. Default: 1.
#'
#' @return A table with the results:
#' \code{VIP} = sqrt(p*q/s), with:
#' \code{p} is the number of variables in each block
#' \code{q} is the explained variance of Y associated to each variable
#' \code{s} is the total Y variance explained by the model
#'
#' @examples
#' MBVIP(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")], 
#'       Y=demo_3_Omics$Y)
#' @export
#' 
MBVIP <- function(data, Y, model = NULL, mc.cores = 1) {
    # Variable format control
    if (!is.list(data)) stop("data is not a list.")
    if (!is.matrix(Y)) stop("Y is not a matrix.")
  
    # Build a model from given data if model is NULL
    if (is.null(model)) {
        model <- RVConsensusOPLS(data=data, Y=Y, modelType="da", A=1)
    }
    if (!is.list(model)) stop("model is not a list.")

    VIP <- parallel::mclapply(1:length(data), mc.cores=mc.cores, function(itable) {
        # Dimensions of the data in the data
        nvariable <- ncol(data[[itable]])
        nsample   <- nrow(model$Model$T)
        ncomp     <- ncol(model$Model$T)

        Qs <- crossprod(Y, model$Model$T) %*% diag(1/diag(crossprod(model$Model$T)), 
                                                   ncol=ncomp) #nclass x ncomp
        Us <- crossprod(t(Y), Qs) %*% diag(1/diag(crossprod(Qs)), 
                                           ncol=ncomp) #nsample x ncomp
        Ws <- crossprod(data[[itable]], Us) %*% diag(1/diag(crossprod(Us)), 
                                                     ncol=ncomp) #nvariable x ncomp
        Ws <- apply(Ws, 2, function(x) x/norm(x, type='2'))
        s <- diag(crossprod(model$Model$T) %*% crossprod(Qs)) # ncomp x ncomp x ncomp x ncomp
    
        VIP.itable <- apply(Ws, 1, function(x) {
            q <- crossprod(s, x^2)
            return (sqrt(nvariable * q / sum(s)))
        })
        return (VIP.itable)
    })

    return (VIP)
}
