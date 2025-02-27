#' @title koplsModel
#' @description Function for training a K-OPLS model. The function constructs a
#' predictive regression model for predicting the values of 'Y' by
#' using the information in 'K'. The explained variation is separated
#' into predictive components (dimensionality determined by the
#' parameter 'A') and 'Y'-orthogonal components (dimensionality determined 
#' by the parameter 'nox').
#'
#' @param K matrix. Kernel matrix (un-centered); K = <phi(Xtr),phi(Xtr)>.
#' @param Y matrix. Response matrix (un-centered/scaled). 
#' @param A numeric. Number of predictive components. Default, 1.
#' @param nox numeric. Number of Y-orthogonal components. Default is 1.
#' @param preProcK character. Pre-processing parameters for the \code{K} matrix:
#' \code{mc} for mean-centering, \code{no} for no centering. Default, 'no'.
#' @param preProcY character. Pre-processing parameters for the \code{Y} matrix:
#' \code{mc} for mean-centering + no scaling, \code{uv} for mc + scaling to unit 
#' variance, \code{pa} for mc + scaling to Pareto, \code{no} for no centering + 
#' no scaling. Default, \code{no}.
#'
#' @returns A list with the following entries:
#' \item{Cp}{ matrix. Y loading matrix.}
#' \item{Sp}{ matrix. The inverse square root of singular values of Y'*K*Y, in
#' diagonal matrix.}
#' \item{Up}{ matrix. Y score matrix.}
#' \item{Tp}{ list. Predictive score matrix for all Y-orthogonal components.}
#' \item{T}{ matrix. Predictive score matrix for the final model.}
#' \item{co}{ list. Y-orthogonal loading vectors.}
#' \item{so}{ list. Eigenvalues from estimation of Y-orthogonal loading vectors.}
#' \item{to}{ list. Weight vector for the i-th latent component of the KOPLS 
#' model.}
#' \item{To}{ matrix. Y-orthogonal score matrix.}
#' \item{toNorm}{ list. Norm of the Y-orthogonal score matrix prior to scaling.}
#' \item{Bt}{ list. T-U regression coefficients for predictions.}
#' \item{A}{ numeric. Number of predictive components.}
#' \item{nox}{ numeric. Number of Y-orthogonal components.}
#' \item{K}{ matrix. The kernel matrix.}
# #' \item{EEprime}{ matrix. The deflated kernel matrix for residual statistics.}
#' \item{sstot_K}{ numeric. Total sums of squares in K.}
#' \item{R2X}{ numeric. Cumulative explained variation for all model components.}
#' \item{R2XO}{ numeric. Cumulative explained variation for Y-orthogonal model 
#' components.}
#' \item{R2XC}{ numeric. Explained variation for predictive model components 
#' after addition of Y-orthogonal model components. Ignored in this version.}
#' \item{sstot_Y}{ numeric. Total sums of squares in Y.}
#' \item{R2Y}{ numeric. Explained variation of Y.}
#' \item{R2Yhat}{ numeric. Variance explained by the i-th latent component of 
#' the model.}
#' \item{preProc$K}{ character. Pre-processing setting for K.}
#' \item{preProc$Y}{ character. Pre-processing setting for Y.}
#' \item{preProc$paramsY}{ character. Pre-processing scaling parameters for Y.}
#'
#' @examples
#' K <- ConsensusOPLS:::koplsKernel(X1 = demo_3_Omics[["MetaboData"]], 
#'                                  X2 = NULL, type='p', params=c(order=1.0))
#' Y <- demo_3_Omics$Y
#' A <- 2
#' nox <- 4
#' preProcK <- "mc"
#' preProcY <- "mc"
#' model <- ConsensusOPLS:::koplsModel(K = K, Y = Y, A = A, nox = nox, 
#'                                     preProcK = preProcK,
#'                                     preProcY = preProcY)
#' ls(model)
#' 
#' @keywords internal
#' @noRd
#' 
koplsModel <- function(K, Y, A = 1, nox = 1,
                       preProcK = "no",
                       preProcY = "no") {
    # Variable format control
    if (!is.matrix(K))
        stop("K is not a matrix.")
    if (!is.matrix(Y))
        stop("Y is not a matrix.")
    if (!is.numeric(A))
        stop("A is not numeric.")
    if (!is.numeric(nox))
        stop("nox is not numeric.")
    preProcK <- match.arg(preProcK,
                          choices = c("mc", "no"),
                          several.ok = F)
    preProcY <- match.arg(preProcY,
                          choices = c("mc", "no", "uv", "pa"),
                          several.ok = F)
    
    # Initialize parameters
    I <- diag(ncol(K))
    
    # Preprocess K
    Kpreproc <- if (preProcK == "mc") koplsCenterKTrTr(KtrTr = K) else K
    # TODO: not necessary to store all K matrices
    Kdeflate <- matrix(data = list(), ncol = nox+1, nrow = nox+1)
    Kdeflate[1,1] <- list(Kpreproc)
    
    # Preprocess Y
    scaleParams <- koplsScale(X = Y, 
                              centerType = ifelse(preProcY == "no",
                                                  "no",
                                                  "mc"),
                              scaleType = ifelse(preProcY == "mc",
                                                 "no",
                                                 preProcY))
    Ypreproc <- scaleParams$X
    
    # KOPLS model estimation
    ## step 1: SVD of Y'KY
    A <- min(A, max(ncol(Y)-1, 1))
    # TOREMOVE: nclass * nsample * nsample * nsample * nsample * nclass = nclass * nclass
    CSV <- svd(x = crossprod(x = Ypreproc, 
                             y = crossprod(x = t(Kdeflate[1,1][[1]]), 
                                           y = Ypreproc)),
               nu = A, nv = A)
    # Extract left singular vectors
    # TOREMOVE: nclass * ncomp, ncomp=A
    Cp <- CSV$u
    rownames(Cp) <- colnames(Ypreproc)
    # Extract the inverse square root of singular values
    Sp <- diag(CSV$d[1:A]^(-1/2), nrow=A)
    
    ## step 2: Define Up
    # TOREMOVE: nsample * nclass * nclass * ncomp = nsample * ncomp
    Up <- crossprod(x = t(Ypreproc), y = Cp)
    
    # Initiate Y-orthogonal related variables
    to <- co <- so <- toNorm <- Tp <- Bt <- list()
    i <- 1
    
    ## step3: Loop over nox iterations
    while (i <= nox) {
        ## step 4: Compute Tp
        # TOREMOVE: nsample * nsample * nsample * ncomp * ncomp * ncomp = nsample * ncomp
        Tp[[i]] <- crossprod(x = Kdeflate[1,i][[1]], 
                             y = tcrossprod(x = Up, y = t(Sp)))
        # TOREMOVE: ncomp * ncomp * ncomp * nsample * nsample * ncomp = ncomp * ncomp
        Bt[[i]] <- crossprod(x = t(solve(crossprod(Tp[[i]]))), 
                             y = crossprod(x = Tp[[i]], y = Up))
        ## step 5: SVD of T'KT
        # TOREMOVE: ncomp * nsample * nsample * nsample * nsample * ncomp = ncomp * ncomp
        temp <- svd(x = crossprod(x = Tp[[i]], 
                                  y = tcrossprod(x = Kdeflate[i,i][[1]] -
                                                     tcrossprod(Tp[[i]]),
                                                 y = t(Tp[[i]]))),
                    nu = 1, nv = 1)
        # TOREMOVE: ncomp * 1
        co[[i]] <- temp$u
        # TOREMOVE: 1
        so[[i]] <- temp$d[1]
        ## TODO: so[[i]] is scalar then?
        
        ## step 6: to
        to[[i]] <- tcrossprod(
            x = tcrossprod(
                x = tcrossprod(
                    x = Kdeflate[i,i][[1]] - tcrossprod(Tp[[i]]), 
                    y = t(Tp[[i]])),
                y = t(co[[i]])), 
            y = t(1/sqrt(so[[i]])))
        
        ## step 7: toNorm
        toNorm[[i]] <- c(sqrt(crossprod(to[[i]])))
        
        ## step 8: Normalize to
        to[[i]] <- to[[i]] / toNorm[[i]]
        
        ## step 9: Update K
        scale_matrix <- I - tcrossprod(to[[i]])
        Kdeflate[1, i+1][[1]] <- tcrossprod(
            x = Kdeflate[1,i][[1]],
            y = t(scale_matrix))
        
        ## step 10: Update Kii
        Kdeflate[i+1, i+1][[1]] <- tcrossprod(
            x = scale_matrix, 
            y = tcrossprod(
                x = t(scale_matrix),
                y = t(Kdeflate[i, i][[1]])))
        
        # Update i
        i <- i + 1
    }## step 11: end loop

    ## step 12: Tp[[nox+1]]
    Tp[[nox+1]] <- crossprod(x = Kdeflate[1, nox+1][[1]], 
                             y = crossprod(x = t(Up), y = Sp))
    colnames(Tp[[nox+1]]) <- paste0("p_", 1:A)
    #TODO: Why is Y-predictive score determined by the optimal number of orthogonal components?
    
    ## step 13: Bt[[nox+1]]
    Bt[[nox+1]] <- crossprod(x = t(solve( crossprod(Tp[[nox+1]]) )),
                             y = crossprod(x = Tp[[nox+1]], y = Up))
    
    # ---------- extra stuff -----------------
    # should work but not fully tested (MB 2007-02-19)
    #sstot_Y <- sum(sum(Ypreproc**2))
    if (max(abs(colMeans(Ypreproc))) > 1e-5)
        stop("Ypreproc is not centered.")
    sstot_Y <- sum(Ypreproc^2)
    #FY <- Ypreproc - tcrossprod(x = Up, y = Cp)
    #R2Y <- 1 - sum(sum(FY^2))/sstot_Y
    # --------- #
    
    #EEprime <- Kdeflate[nox+1, nox+1][[1]] - tcrossprod(Tp[[nox+1]])
    sstot_K <- sum(diag(Kdeflate[1,1][[1]]))
    
    # Explained variance
    R2X <- 1 - sapply(
        1:(nox+1),
        function(i) {
            sum(diag(Kdeflate[i,i][[1]] - tcrossprod(Tp[[nox+1]])))
        })/sstot_K
    R2XC <- rep(
        1 - sum(
            diag(Kdeflate[1,1][[1]] - tcrossprod(Tp[[nox+1]]))
        )/sstot_K,
        times = nox+1)
    R2XO <- 1 - sapply(1:(nox+1),
                       function(i) {
                           sum(diag(Kdeflate[i,i][[1]]))    
                       })/sstot_K
    Yhat <- lapply(1:(nox+1),
                   function(i) {
                       crossprod(t(Tp[[i]]), tcrossprod(Bt[[i]], Cp))
                   })
    R2Yhat <- 1 - sapply(Yhat,
                         function(yh) {
                             #sum(sum((yh - Ypreproc)^2))/sstot_Y 
                             sum((yh - Ypreproc)^2)/sstot_Y
                         })
    names(R2X) <- names(R2XC) <- names(R2XO) <- names(R2Yhat) <-
        c("p", paste0("po", 1:nox))
    
    # Convert to matrix structure
    if (nox > 0) {
        To <- do.call(cbind, to)
        colnames(To) <- paste0("o_", 1:ncol(To))
    } else {
        To <- NULL
    }
    
    # Group parameters
    params <- list("nPcomp"   = A,
                   "nOcomp"   = nox, 
                   "sstot_K"  = sstot_K,
                   "sstot_Y"  = sstot_Y,
                   "preProcK" = preProcK, 
                   "preProcY" = preProcY)
    
    koplsObj <- list("params"   = params,
                     "scoresP"  = as.matrix(Tp[[nox+1]]),
                     "scoresO"  = as.matrix(To),
                     "Cp"       = Cp,
                     "Sp"       = Sp,
                     "Up"       = Up,
                     "Tp"       = Tp,
                     "co"       = co,
                     "so"       = so,
                     "to"       = to,
                     "toNorm"   = toNorm,
                     "Bt"       = Bt,
                     "K"        = Kdeflate,
                     
                     #extra stuff
                     #"EEprime"  = EEprime,
                     "R2X"      = R2X,
                     "R2XO"     = R2XO,
                     #"R2Y"     = R2Y,
                     "R2Yhat"   = R2Yhat, # R2Yhat 22 Jan 2010 / MR
                     
                     #pre-processing
                     "preProc"  = list("paramsY" = scaleParams)
    )
    class(koplsObj) <- c("kopls")
    
    return (koplsObj)
}
