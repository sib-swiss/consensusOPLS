#' @title koplsCrossValSet
#' @description Generates set of training/test observations for the 
#' Cross-Validation.
#'
#' @param K matrix. Kernel matrix. 
#' @param Y matrix. Response matrix.
#' @param cvFrac numeric. Fraction (in percent) of observations used in the 
#' training set. Default, 4/5.
#' @param cvType character. Type of cross-validation wanted. It can be \code{nfold} 
#' for n-fold, \code{mccv} for Monte Carlo CV, \code{mccvb} for Monte Carlo 
#' class-balanced CV. Default is \code{nfold}.
#' @param nfold numeric. Number of total nfold rounds (if cvType = \code{nfold}).
#' @param nfoldRound numeric. Current nfold rounds (if cvType = \code{nfold}).
#' @param mc.cores Number of cores for parallel computing. Default: 1.
#'
#' @returns A list with the following entries:
#' \item{CV_param}{ data frame. It contains \code{cvType} a character for the 
#' Cross-validation type, \code{nfold} the total number of nfold, 
#' \code{nfoldRound} the current nfold round and \code{class} the object class
#' wich is \code{crossValSet}.}
#' \item{KTrTr}{ matrix. Kernel training matrix; KTrTr = <phi(Xtr),phi(Xtr)>.}
#' \item{KTeTr}{ matrix. Kernel test/training matrix; KTeTr = <phi(Xte),phi(Xtr)>.}
#' \item{KTeTe}{ matrix. Kernel test matrix; KTeTe = <phi(Xte),phi(Xte)>.}
#' \item{yTraining}{ matrix. Y training set.}
#' \item{yTest}{ matrix. Y test set.}
#' \item{trainingIndex}{ integer. Indices of training set observations.}
#' \item{testIndex}{ numeric. Indices of test set observations.}
#'
#' @examples
#' Y <- matrix(stats::rnorm(n = 28), nrow = 14)
#' K <- matrix(stats::rnorm(n = 140), nrow = 14)
#' cvType <- "nfold"  
#' cvFrac <- 0.75 
#' nfold <- 5  
#' nfoldRound <- 1 
#' cvs <- ConsensusOPLS:::koplsCrossValSet(K = K, Y = Y, 
#'                                          cvFrac = cvFrac,
#'                                          cvType = cvType, nfold = nfold, 
#'                                          nfoldRound = nfoldRound)
#' cvs 
#' 
#' @import parallel
#' @keywords internal
#' @noRd
#' 
koplsCrossValSet <- function(K, Y, cvFrac = 4/5, cvType = "nfold", 
                             nfold = NA, nfoldRound = NA, mc.cores = 1, 
                             random.seed = 10403) {
    # Variable format control
    if (!is.matrix(K)) stop("K is not a matrix.")
    if (!is.matrix(Y)) stop("Y is not a matrix.")
    if (!is.numeric(cvFrac)) stop("cvFrac is not numeric.")
    if (!is.character(cvType)) stop("cvType is not a character.")
    else if (!(cvType %in% c("nfold", "mccv", "mccvb")))
        stop("cvType must be `nfold`, `mccv` or `mccvb`.")
    if (!is.na(nfold)) {
        if (!is.numeric(nfold)) stop("nfold is not a number.")
    } else if (cvType != "nfold")
        warning("cvType is not nfold, nfold is not defined (missing argument).")
    
    if (!is.na(nfoldRound)) {
        if (!is.numeric(nfoldRound)) stop("nfoldRound is not a number.")
    } else if (cvType != "nfold")
        warning("cvType is not nfold, nfoldRound is not defined (missing argument).")
    
    # Set random seed
    #set.seed(random.seed)
    
    # Define Monte-Carlos Cross Validation - class Balanced
    if (cvType == "mccvb") {
        # check if Y is dummy or labels
        if (all(unique(x = Y) %in% c(0, 1))) {
            classVect <- koplsReDummy(Y = Y)
        } else{
            classVect <- Y
        }
        
        # Find all class labels
        uniqueClass <- unique(x = classVect)
        
        # Find samples of each class
        indList <- mclapply(X = uniqueClass, 
                            mc.cores = mc.cores,
                            FUN = function(i) {
                                ind <- which(classVect == i)
                                rand_ind <- sample(x = ind)
                                return (sample(x = ind))
                            })
        
        # Combine indices for all classes
        if (nfoldRound > 0) {
            trainInd <- unlist(x = indList)
            predInd <- setdiff(x = seq_len(nrow(Y)), y = trainInd)
        } else {
            trainInd <- predInd <- 1:nrow(Y)
        }
    }
    
    # Define Monte-Carlos Cross Validation
    if (cvType == "mccv") {
        # Divides the sample randomly into train and test
        if (nfoldRound > 0) {
            trainInd <- sample.int(nrow(Y), floor(x = nrow(Y)*cvFrac))
            predInd <- setdiff(1:nrow(Y), trainInd)
        } else {
            trainInd <- predInd <- 1:nrow(Y)
        }
    }

    # Define N-fold Cross Validation
    if (cvType == "nfold") {
        if (nfoldRound > 0) {
            predInd <- seq.int(from=nfoldRound, to=nrow(Y), by=nfold)
            trainInd <- setdiff(1:nrow(Y), predInd)
        } else {
            trainInd <- predInd <- 1:nrow(Y)
        }
    }
    
    # Construct Kernel/Y matrices for training/test
    if (nrow(K) == ncol(K)) {
        KTrTr <- K[trainInd, trainInd, drop=FALSE]
        KTeTr <- K[predInd, trainInd, drop=FALSE]
        KTeTe <- K[predInd, predInd, drop=FALSE]
    } else {
        KTrTr <- NA
        KTeTr <- NA
        KTeTe <- NA
    }
    
    # Group parameters in data.frame
    CV_param <- data.frame("cvType" = cvType,
                           "nfold" = nfold,
                           "nfoldRound" = nfoldRound,
                           "class" = "crossValSet")
    
    # Return the final CV Set
    return (list("CV_param" = CV_param,
                 "KTrTr" = KTrTr,
                 "KTeTr" = KTeTr,
                 "KTeTe" = KTeTe,
                 "yTraining" = Y[trainInd, , drop=FALSE],
                 "yTest" = Y[predInd, , drop=FALSE],
                 "trainingIndex" = trainInd,
                 "testIndex" = predInd))
}
