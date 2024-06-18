#' @title Variable Importance in Projection
#' @description Calculate the VIP (Variable Importance in Projection) for each
#' variable in a \code{ConsensusOPLS} model.
#' 
#' @param data A list of data blocks.
#' @param Y A vector, factor, dummy matrix or numeric matrix for the response.
#' @param model An object inheriting from class \code{ConsensusOPLS}, 
#' representing a Consensus OPLS fitted model. Default, NULL, a model will be
#' constructed.
#' @param combine A logical indicating the formula of VIP. Default, FALSE, the
#' formula 2 in Galindo-Prieto et al. (2013 )Variable influence on projection for
#' orthogonal projections to latent structures is used. If TRUE, formula 4.
#' @param ... \code{RVConsensusOPLS} arguments.
#'
#' @returns A list of VIP for variables in the data blocks:
#' \code{VIP = sqrt(p*q/s)}, where
#' \code{p} is the number of variables in each block,
#' \code{q} the explained variance of Y associated to each variable, and 
#' \code{s} the total Y variance explained by the model.
#'
#' @examples
#' vip <- ConsensusOPLS:::VIP(
#'     data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")], 
#'     Y=demo_3_Omics$Y)
#' str(vip)
#' @keywords internal
#' @noRd
#' 
VIP <- function(data, Y, model = NULL, combine = F, ...) {
    # Variable format control
    if (!is.list(data)) stop("data is not a list.")
    if (!is.matrix(Y)) stop("Y is not a matrix.")
    
    # Build a model from given data if model is NULL
    if (is.null(model)) {
        rvcopls <- RVConsensusOPLS(data=data, Y=Y, ...)
        model <- rvcopls$Model
    }
    if (!is.list(model)) stop("model is not a list.")
    
    if (combine) {
        stop("To explore")
    } else {
        # compute VIP for predictive and orthogonal parts
        VIPs <- lapply(c("scoresP", "scoresO"), function(sc) {
            VIP.sc <- lapply(1:length(data), function(itable) {
                # Dimensions of the data in the data
                nvariable <- ncol(data[[itable]])
                nsample   <- nrow(model[[sc]])
                ncomp     <- ncol(model[[sc]])
                
                # nclass x ncomp
                Qs <- crossprod(Y, model[[sc]]) %*% diag(1/diag(crossprod(model[[sc]])), 
                                                          ncol=ncomp)
                # nsample x ncomp
                Us <- crossprod(t(Y), Qs) %*% diag(1/diag(crossprod(Qs)), 
                                                   ncol=ncomp)
                # nvariable x ncomp
                # NB. why is it not model$loadings?
                Ws <- crossprod(data[[itable]], Us) %*% diag(1/diag(crossprod(Us)), 
                                                             ncol=ncomp)
                Ws <- apply(Ws, 2, function(x) x/norm(x, type='2'))
                # ncomp x ncomp x ncomp x ncomp
                s <- diag(crossprod(model[[sc]]) %*% crossprod(Qs))
                
                VIP.itable <- do.call(rbind, lapply(1:nrow(Ws), function(j) {
                    q <- crossprod(s, Ws[j,]^2)
                    return (sqrt(nvariable * q / sum(s)))
                }))
                rownames(VIP.itable) <- rownames(Ws)
                
                return (VIP.itable)
            })
            names(VIP.sc) <- names(data)
            return (VIP.sc)
        })
        names(VIPs) <- c("scoresP", "scoresO")
        
        # reform the VIP of each data block
        VIP <- lapply(1:length(data), function(itable) {
            VIP.itable <- data.frame(p=VIPs$scoresP[[itable]], 
                                     o=VIPs$scoresO[[itable]])
            VIP.itable$tot <- sqrt((VIP.itable$p^2 + VIP.itable$o^2)/2)
            
            return (VIP.itable)
        }) 
        names(VIP) <- names(data)
    }
    
    return (VIP)
}
