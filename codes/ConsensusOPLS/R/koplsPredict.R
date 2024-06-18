#' @title koplsPredict
#' @description Performs prediction of new samples from an existing K-OPLS model.
#' The function projects the Y-predictive and Y-orthogonal scores components 
#' to predict a value of the response matrix Y. The dimensions of the 
#' parameters is determined from the specified model.
#'
#' @param KteTr matrix. The hybrid test/training kernel matrix. 
#' KteTr = <phi(Xte),phi(Xtr)>.
#' @param Ktest matrix. The pure test kernel matrix. Ktest = <phi(Xte),phi(Xte)>.
#' @param Ktrain matrix. The training kernel matrix (same as used in model 
#' training). Ktrain = <phi(Xtr),phi(Xtr)>.
#' @param model list. K-OPLS model object.
#' @param nox numeric. Number of Y-orthogonal components. If not specified, the 
#' number used during model training will be employed.
#' @param rescaleY logical. If \code{TRUE}, predicted values of the
#' response (Yhat) is rescaled according to the pre-processing settings of 
#' the model. If \code{FALSE}, Yhat is not rescaled (default).
#'
#' #TODO: Why Tp, to, T, ... are produced here?
#' @returns A list with the following entries:
#' \item{Tp}{ matrix. Predicted predictive score matrix for all generations 0: 
#' \code{nox} of Y-orthogonal vectors.}
#' \item{to}{ vector. Predicted Y-orthogonal score vectors.}
#' \item{T}{ matrix. Predictive score matrix for the final model.}
#' \item{KteTr}{ matrix. Predictive score matrix for the final model with 
#' \code{nox} Y-orthogonal vectors.}
#' \item{EEprime}{ matrix. Calculated residuals for the test kernel \code{Ktest}, 
#' useful e.g. for residual statistics.}
#' \item{Yhat}{ matrix. Predicted values of the response matrix.}
#'
#' @examples
#' Xte <- matrix(data = stats::rnorm(n = 20), ncol=5)
#' Xtr <- matrix(data = stats::rnorm(n = 25), ncol=5)
#' KteTe <- ConsensusOPLS:::koplsKernel(X1 = Xte, X2 = Xte, 
#'                                      type='p', params=c(order=1.0))
#' KteTr <- ConsensusOPLS:::koplsKernel(X1 = Xte, X2 = Xtr, 
#'                                      type='p', params=c(order=1.0))
#' KtrTr <- ConsensusOPLS:::koplsKernel(X1 = Xtr, X2 = Xtr, 
#'                                      type='p', params=c(order=1.0))
#' 
#' Y <- matrix(data = stats::rnorm(n = 5), nrow = 5)
#' A <- 2
#' nox <- 4
#' preProcK <- "mc"
#' preProcY <- "mc"
#' model <- ConsensusOPLS:::koplsModel(K = KtrTr, Y = Y, A = A, nox = nox, 
#'                                     preProcK = preProcK, preProcY = preProcY)
#' pred <- ConsensusOPLS:::koplsPredict(KteTr = KteTr, Ktest = KteTe, 
#'                                      Ktrain = KtrTr, model = model, nox = nox,
#'                                      rescaleY = FALSE)
#' pred
#' @keywords internal
#' @noRd
#' 
koplsPredict <- function(KteTr, Ktest, Ktrain,
                         model, nox, rescaleY = FALSE) {
    # Variable format control
    if (!is.matrix(KteTr) || !is.matrix(Ktest) || !is.matrix(Ktrain)) {
        stop("One or more kernel inputs are not matrices.")
    }
    if (!is.list(model)) stop("model is not a list containing model parameters.")
    else if (model$params$class != "kopls") stop("Model must be of type `kopls`.")
    if (!is.null(nox)) {
        if (!is.numeric(nox)) stop("nox is not numeric.")
        if (nox > model$params$nOcomp) {
            warning("Number of Y-orthogonal components to use is higher than in model.
              Setting number of Y-orthogonal to max in model.")
            nox <- model$params$nOcomp
        }
    } else stop('Number of Y-orthogonal components to use is missing.')
    if (is.null(rescaleY)) rescaleY <- 0
    else if (!is.logical(rescaleY)) stop("rescaleY is not logical.")
    
    # Step1: mean centering of K matrices
    # the order of the code below is important
    KteTepreproc <- if (model$params$preProcK == "mc") 
        koplsCenterKTeTe(KteTe = Ktest, KteTr = KteTr, KtrTr = Ktrain) else 
            Ktest

    KteTedeflate <- matrix(data = list(NULL), 
                           nrow = model$params$nOcomp + 1, 
                           ncol = model$params$nOcomp + 1)
    KteTedeflate[1,1][[1]] <- KteTepreproc
    
    KteTrpreproc <- if (model$params$preProcK == "mc") 
        koplsCenterKTeTr(KteTr = KteTr, KtrTr = Ktrain) else 
            KteTr
    KteTrdeflate <- matrix(data = list(NULL), 
                           nrow = model$params$nOcomp + 1, 
                           ncol = model$params$nOcomp + 1)
    KteTrdeflate[1,1][[1]] <- KteTrpreproc

    # Initialize variables
    to <- Tp <- list()
    
    # Step2: KOPLS prediction
    ## Step2.1: for each Y-orth component
    if (nox > 0) {
        for (i in 1:nox) {
            ## Step2.2: Predicted predictive score matrix
            Tp[[i]] <- crossprod(x = t(KteTrdeflate[i,1][[1]]),
                                 y = tcrossprod(x = model$Up, 
                                                y = t(model$Sps)))
            
            # Step2.3: Predicted Y-orthogonal score vectors
            # to[[i]] <- crossprod(x = t((KteTrdeflate[i,i][[1]] - 
            #                                 tcrossprod(x = Tp[[i]], 
            #                                            y = model$Tp[[i]]))),
            #                      y = tcrossprod(x = model$Tp[[i]], 
            #                                     y = tcrossprod(x = t(sqrt(model$so[[i]])),
            #                                                    #y = t(model$co[[i]]))))
            #                                                    y = model$co[[i]])))
            to[[i]] <- tcrossprod(x = tcrossprod(x = tcrossprod(x = KteTrdeflate[i,i][[1]] - tcrossprod(x = Tp[[i]], 
                                                                                                    y = model$Tp[[i]]), 
                                                                y = t(model$Tp[[i]])),
                                                 y = t(model$co[[i]])),
                                  y = t(1/sqrt(model$so[[i]])))
            
            # Step2.4: Normalize to
            to[[i]] <- to[[i]]/model$toNorm[[i]]
            
            # Step2.4.5: deflate KteTedeflate (this is an EXTRA feature - not in alg. in paper)
            KteTedeflate[i+1,i+1][[1]] <- KteTedeflate[i,i][[1]] - 
                crossprod(x = t(KteTrdeflate[i,i][[1]]), 
                          y = tcrossprod(x = model$to[[i]], 
                                         y = to[[i]])) - 
                crossprod(x = t(to[[i]]), 
                          y = crossprod(x = model$to[[i]], 
                                        y = t(KteTrdeflate[i,i][[1]]))) +
                crossprod(x = t(to[[i]]),
                          y = crossprod(x = model$to[[i]],
                                        y = crossprod(x = t(model$K[i,i][[1]]),
                                                      y = tcrossprod(x = model$to[[i]],
                                                                     y = to[[i]]))))
            
            # Step2.5: Update KteTrdeflate
            KteTrdeflate[i+1,1][[1]] <- KteTrdeflate[i,1][[1]] - 
                crossprod(x = t(to[[i]]),
                          y = crossprod(x = model$to[[i]], 
                                        y = t(model$K[1,i][[1]])))
            
            # Step2.6: Update KteTrdeflate
            KteTrdeflate[i+1,i+1][[1]] <- KteTrdeflate[i,i][[1]] - 
                crossprod(x = t(KteTrdeflate[i,i][[1]]),
                          y = tcrossprod(model$to[[i]])) - 
                crossprod(x = t(to[[i]]),
                          y = crossprod(x = model$to[[i]], 
                                        y = model$K[i,i][[1]])) + 
                crossprod(x = t(to[[i]]), 
                          y = crossprod(x = model$to[[i]], 
                                        y = crossprod(x = t(model$K[i,i][[1]]), 
                                                      y = tcrossprod(model$to[[i]]))))
        } # Step2.7: end loop
    } else if (nox == 0) {
        i <- 0
    }
    
    Tp[[i+1]] <- crossprod(x = t(KteTrdeflate[i+1,1][[1]]),
                           y = tcrossprod(model$Up, t(model$Sps)))
    # TOREMOVE: nsample * ncomp * ncomp * ncomp * ncomp * nclass = nsample * nclass
    Yhat <- crossprod(x = t(Tp[[i+1]]),
                      y = tcrossprod(model$Bt[[i+1]], model$Cp))
    
    if (!is.null(rescaleY)) {
        if (model$params$preProcY == "no") {
            if (length(model$preProc$paramsY) == 1 && model$preProc$paramsY == "no") {
                scaleParams <- list()
                scaleParams$centerType <- "no"
                scaleParams$scaleType <- "no"
            } else {
                scaleParams <- model$preProc$paramsY
            }
            YhatRescaled <- koplsRescale(scaleS = scaleParams, varargin = Yhat)
            Yhat <- YhatRescaled$X
        } else {
            warning("Attempted re-scale of Yhat although no pre-processing 
              parameters have been set.")
        }
    }
    
    #---- Extra stuff ----------------------------------
    EEprime <- KteTedeflate[i+1,i+1][[1]] - tcrossprod(Tp[[i+1]])
    #--------------------------------------------------
    
    # Return the list of prediction parameters
    return (list(#"Tp"           = Tp,
                 #"to"           = to,
                 #"T"            = Tp[[nox+1]],
                 #"KteTrdeflate" = KteTrdeflate,
                 #"EEprime"      = EEprime,
                 "Yhat"         = Yhat))
}
