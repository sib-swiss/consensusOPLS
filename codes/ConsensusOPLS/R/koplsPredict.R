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
#' @return A list with the following entries:
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
#' Xte <- matrix(rnorm(20), ncol=5)
#' Xtr <- matrix(rnorm(25), ncol=5)
#' KteTe <- ConsensusOPLS:::koplsKernel(Xte, Xte, Ktype='g', params=c(sigma=1.0))
#' KteTr <- ConsensusOPLS:::koplsKernel(Xte, Xtr, Ktype='g', params=c(sigma=1.0))
#' KtrTr <- ConsensusOPLS:::koplsKernel(Xtr, Xtr, Ktype='g', params=c(sigma=1.0))
#' 
#' Y <- matrix(rnorm(15), nrow = 5)
#' A <- 2
#' nox <- 4
#' preProcK <- "mc"
#' preProcY <- "mc"
#' model <- ConsensusOPLS:::koplsModel(K = KtrTr, Y = Y, A = A, nox = nox, 
#'                                     preProcK = preProcK, preProcY = preProcY)
#' ##pred <- ConsensusOPLS:::koplsPredict(KteTr, KteTe, KtrTr, model, nox)
#'              

#' @keywords internal

koplsPredict <- function(KteTr, Ktest, Ktrain,
                         model, nox, rescaleY = FALSE) {
    # Variable format control
    if (!is.matrix(KteTr) || !is.matrix(Ktest) || !is.matrix(Ktrain)) {
        stop("One or more kernel inputs are not matrices.")
    }
    if (!is.list(model)) stop("model is not a list containing model parameters.")
    else if(model$class != "kopls") stop("Model must be of type `kopls`.")
    if(!is.null(nox)) {
        if (!is.numeric(nox)) stop("nox is not numeric.")
        if (nox > model$nox) {
            warning("Number of Y-orthogonal components to use is higher than in model.
              Setting number of Yorth to max in model.")
            nox <- model$nox
        }
    } else {
        stop('Number of Y-orthogonal components to use is missing.')
    }
    if (is.null(rescaleY)) rescaleY <- 0
    else if (!is.logical(rescaleY)) stop("rescaleY is not logical.")
    
    # Step1: mean centering of K matrices
    # the order of the code below is important
    KteTeMc <- Ktest
    if (model$preProc$K == "mc") {
        KteTeMc <- koplsCenterKTeTe(KteTe = Ktest, 
                                    KteTr = KteTr, 
                                    KtrTr = Ktrain)
    }
    KteTe <- matrix(data = list(NULL), nrow = model$nox+1, ncol = model$nox+1)
    KteTe[1,1][[1]] <- KteTeMc
    
    KteTrMc <- KteTr
    if (model$preProc$K == "mc") {
        KteTrMc <- koplsCenterKTeTr(KteTr = KteTr, KtrTr = Ktrain)
    }
    KteTrTmp <- matrix(data = list(NULL), nrow = model$nox+1, 
                       ncol = model$nox+1)
    KteTrTmp[1,1][[1]] <- KteTrMc
    KteTr <- KteTrTmp
    
    # Initialize variables
    to <- list()
    Tp <- list()
    
    # Step2: KOPLS prediction
    ## Step2.1: for each Y-orth component
    if (nox > 0) {
        for (i in 1:nox) {
            ## Step2.2: Predicted predictive score matrix
            Tp[[i]] <- tcrossprod(KteTr[i,1][[1]],
                                  tcrossprod( t(model$Sps), t(model$Up))) # TODO:check dimensions
            
            # Step2.3: Predicted Y-orthogonal score vectors
            to[[i]] <- crossprod(t((KteTr[i,i][[1]] - tcrossprod(Tp[[i]], 
                                                                 model$Tp[[i]]))),
                                 tcrossprod(model$Tp[[i]], 
                                            tcrossprod(t(sqrt(model$so[[i]])), 
                                                       t(model$co[[i]]))))
            
            # Step2.4: Normalize to
            to[[i]] <- to[[i]]/model$toNorm[[i]]
            
            # Step2.4.5: deflate KteTe (this is an EXTRA feature - not in alg. in paper)
            KteTe[i+1,i+1][[1]] <- KteTe[i,i][[1]] - 
                crossprod(t(KteTr[i,i][[1]]), tcrossprod(model$to[[i]], to[[i]])) - 
                crossprod(t(to[[i]]), crossprod(model$to[[i]], t(KteTr[i,i][[1]]))) +
                crossprod(t(to[[i]]),
                          crossprod(model$to[[i]],
                                    crossprod(t(model$K[i,i][[1]]),
                                              tcrossprod(model$to[[i]],
                                                         to[[i]]))))
            
            # Step2.5: Update KTeTr
            KteTr[i+1,1][[1]] <- KteTr[i,1][[1]] - 
                crossprod(t(to[[i]]),
                          crossprod(model$to[[i]], t(model$K[1,i][[1]])))
            
            # Step2.6: Update KTeTr
            KteTr[i+1,i+1][[1]] <- KteTr[i,i][[1]] - 
                crossprod(t(KteTr[i,i][[1]]),
                          tcrossprod(model$to[[i]])) - 
                crossprod(t(to[[i]]),
                          crossprod(model$to[[i]], model$K[i,i][[1]])) + 
                crossprod(t(to[[i]]), 
                          crossprod(model$to[[i]], 
                                    crossprod(t(model$K[i,i][[1]]), 
                                              tcrossprod(model$to[[i]]))))
        } # Step2.7: end loop
    }
    
    if (nox == 0) i <- 0
    
    Tp[[i+1]] <- tcrossprod(KteTr[i+1,1][[1]],
                            tcrossprod(t(model$Sps), t(model$Up)))
    Yhat <- crossprod(t(Tp[[i+1]]),
                      tcrossprod(model$Bt[[i+1]], model$Cp))
    
    if (!is.null(rescaleY)) {
        if (model$preProc$Y == "no") {
            if(length(model$preProc$paramsY) == 1 & model$preProc$paramsY == "no") {
                scaleParams <- list()
                scaleParams$centerType <- "no"
                scaleParams$scaleType <- "no"
            } else {
                scaleParams <- model$preProc$paramsY
            }
            YhatRescaled <- koplsRescale(scaleS = scaleParams, 
                                         varargin = Yhat)
            Yhat <- YhatRescaled$X
        } else {
            warning("Attempted re-scale of Yhat although no pre-processing 
              parameters have been set.")
        }
    }
    
    #---- Extra stuff ----------------------------------
    EEprime <- KteTe[i+1,i+1][[1]] - tcrossprod(Tp[[i+1]])
    #--------------------------------------------------
    
    # Return the list of prediction parameters
    return (list("Tp" = Tp,
                 "to" = to,
                 "T" = Tp[[nox+1]],
                 "KteTr" = KteTr,
                 "EEprime" = EEprime,
                 "Yhat" = Yhat))
}
