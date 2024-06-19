#' @title ConsensusOPLSCV
#' @description Function for performing Kernel-OPLS cross-validation for a set 
#' of Y-orthogonal components.
#' This function was first implemented in KOPLS1.1 MATLAB package, but was
#' modified in 2012 for the ConsensusOPLS method.
#'
#' @param K matrix. The kernel matrix (un-centered); see \code{koplsKernel} for 
#' details.
#' @param Y matrix. The response matrix (un-centered/scaled). Could be binary 
#' (for discriminant analysis) or real-valued (for classical OPLS analysis).
#' @param maxPcomp numeric. The number of Y-predictive components. 
#' @param maxOcomp numeric. The number of Y-orthogonal components.
#' @param nfold numeric. Number of cross-validation rounds (integer).
#' @param cvType character. Type of cross-validation used. Either \code{nfold} 
#' for n-fold cross-validation, \code{mccv} for Monte Carlo cross-validation or 
#' \code{mccvb} for Monte Carlo class-balanced cross-validation. See also 
#' \code{koplsCrossValSet} for details. 
#' @param preProcK character. Pre-processing settings for the kernel matrix. 
#' Either \code{mc} for mean-centering or \code{no} for no pre-processing.
#' @param preProcY character. Pre-processing parameter for Y. Either \code{mc} 
#' for mean-centering, \code{uv} for mean-centering and scaling to unit-variance, 
#' \code{pareto} for mean-centering and Pareto-scaling or \code{no} for no 
#' mean-centering and no scaling.
#' @param cvFrac numeric. Fraction of observations in the training set during 
#' cross-validation. 
#' @param modelType character. Type of model used for the ConsensusOPLS method. 
#' It could be \code{da} for discriminant analysis or \code{reg} for regression. 
#' If \code{da} (default), sensitivity and specificity will be calculated.
#' @param verbose logical. Indicates whether the user wants to see the 
#' progress bar printed. If \code{FALSE}, no output will be printed. If 
#' \code{TRUE} (default) some output will be printed regarding the 
#' cross-validation progress.
#'
#' @returns A list with the following entries (diagnostic parameters which can be 
#' used to determine the optimal number of model components):
#' \item{Model}{ list. The training a K-OPLS model. It contrains:}
#'      \item{Cp}{ matrix. Y loading matrix.}
#'      \item{Sp}{ matrix. Sigma matrix, containing singular values from 
#'      \code{t(Y)* K *Y} used for scaling.}
#'      \item{Sps}{ matrix. Scaled Sigma matrix, containing scaled singular 
#'      values.}
#'      \item{Up}{ matrix. Y score matrix.}
#'      \item{Tp}{ list. Predictive score matrix for all Y-orthogonal components.}
#'      \item{T}{ matrix. Predictive score matrix for the final model.}
#'      \item{co}{ list. Y-orthogonal loading vectors.}
#'      \item{so}{ list. Eigenvalues from estimation of Y-orthogonal loading 
#'      vectors.}
#'      \item{to}{ list. Weight vector for the i-th latent component of the 
#'      KOPLS model.}
#'      \item{To}{ matrix. Y-orthogonal score matrix.}
#'      \item{toNorm}{ list. Norm of the Y-orthogonal score matrix prior to 
#'      scaling.}
#'      \item{Bt}{ list. T-U regression coefficients for predictions.}
#'      \item{A}{ numeric. Number of predictive components.}
#'      \item{nox}{ numeric. Number of Y-orthogonal components.}
#'      \item{K}{ matrix. The kernel matrix.}
#'      \item{EEprime}{ matrix. The deflated kernel matrix for residual 
#'      statistics.}
#'      \item{sstot_K}{ numeric. Total sums of squares in \code{K}.}
#'      \item{R2X}{ numeric. Cumulative explained variation for all model 
#'      components.}
#'      \item{R2XO}{ numeric. Cumulative explained variation for Y-orthogonal 
#'      model components.}
#'      \item{R2XC}{ numeric. Explained variation for predictive model components 
#'      after addition of Y-orthogonal model components.}
#'      \item{sstot_Y}{ numeric. Total sums of squares in Y.}
#'      \item{R2Y}{ numeric. Explained variation of Y.}
#'      \item{R2Yhat}{ numeric. Variance explained by the i-th latent component 
#'      of the model.}
#'      \item{preProc$K}{ character. Pre-processing setting for K.}
#'      \item{preProc$Y}{ character. Pre-processing setting for Y.}
#'      \item{preProc$paramsY}{ character. Pre-processing scaling parameters for 
#'      Y.}
#' \item{cv}{ list. The cross-validation results.}
#'      \item{Yhat}{ matrix. Predicted Y values.}
#'      \item{AllYhat}{ matrix. All predicted Y values as a concatenated matrix.}
#'      \item{Tcv}{ matrix. Predictive score vector T for all cross-validation 
#'      rounds.}
#'      \item{Q2Yhat}{ matrix. Total Q-square result for all Y-orthogonal 
#'      components.}
#'      \item{Q2YhatVars}{ matrix. Q-square result per Y-variable for all 
#'      Y-orthogonal components.}
#'      \item{cvTestIndex}{ matrix. Indices for the test set observations during 
#'      the cross-validation rounds.}
#'      \item{cvTrainIndex}{ matrix. Indices for the training set observations 
#'      during the cross-validation rounds.}
#' \item{da}{ list. Cross-validation results specifically for discriminant 
#' analysis:}
#'      \item{totalResults}{ list. The results over all classes and averages.}
#'      \item{classResults}{ list. The results for each class.}
#'      \item{predClass}{ integer. Predicted class list per class and Y-orthogonal
#'      components.}
#'      \item{trueClass}{ integer. Predicted class list per class and Y-orthogonal 
#'      components.}
#'      \item{sensSpec}{ integer. Sensitivity and specificity values per class and
#'      Y-orthogonal components.}
#'      \item{confusionMatrix}{ matrix. Confusion matrix during cross-validation
#'      rounds.}
#'      \item{nclasses}{ integer. Number of classes in model.}
#'      \item{decisionRule}{ character. Decision rule used: \code{max} or 
#'      \code{fixed}.}
#'      \item{args}{ list. Arguments to the function:}
#'            \item{maxOcomp}{ integer. Number of Y-orthogonal components.}
#'            \item{maxPcomp}{ integer. Number of Y-predictive components.}
#' \item{class}{ character. Model class is \code{koplscv}.}
#' @import utils
#' @keywords internal 
#' @noRd
#' 
ConsensusOPLSCV <- function(K, Y, maxPcomp, maxOcomp, 
                            modelType = "da", 
                            cvType = 'nfold',
                            nfold = 5,
                            nMC = 100,
                            cvFrac = 4/5,
                            preProcK = "no",
                            preProcY = "no",
                            mc.cores = 1, 
                            verbose = FALSE) {
    # ----- Variable format control (part 1)
    if (!is.matrix(K)) stop("K is not a matrix.")
    if (!is.numeric(maxPcomp)) stop("maxPcomp is not numeric.")
    if (!is.numeric(maxOcomp)) stop("maxOcomp is not numeric.")
    if (!is.numeric(nfold)) stop("nfold is not numeric.")
    if (!is.numeric(nMC)) stop("nfold is not numeric.")
    if (!is.character(cvType)) stop("cvType is not a character.")
    else if(!(cvType %in% c("nfold", "mccv", "mccvb")))
        stop("cvType must be `nfold`, `mccv` or `mccvb`.")
    if (!is.character(preProcK)) stop("preProcK is not a character.")
    else if(!(preProcK %in% c("mc", "no"))) stop("preProcK must be `mc` or `no`.")
    if (!is.character(preProcY)) stop("preProcY is not a character.")
    else if(!(preProcY %in% c("mc", "uv", "pareto", "no")))
        stop("preProcY must be `mc`, `uv`, `pareto` or `no`.")
    if (!is.numeric(cvFrac)) stop("cvFrac is not numeric.")
    if (!is.character(modelType)) stop("modelType is not a character.")
    else if(!(modelType %in% c("da", "reg"))) stop("modelType must be `da` or `reg`.")
    if (!is.logical(verbose)) stop("verbose must be `TRUE` or `FALSE`.")
    
    # ----- Default values control
    if (cvType == "mccvb" && modelType != "da")
        stop("Class balanced MC CV only applicable to `da` modelling.")

    # ----- Variable format control (part 2)
    if (modelType == "da") {
        # Check the response matrix
        if (! (all(Y %in% c(0, 1)) && ncol(Y) > 1)) stop("Y is not a dummy matrix.")
        
        classVect <- koplsReDummy(Y = Y)
        nclasses <- length(unique(x = classVect))
    } else {
        if (! (is.numeric(Y) && (ncol(Y) == 1 || is.vector(Y)))) {
            stop("Y should be a numeric vector.")
        }
    }

    # ----- Convert Y-scaling to more explicit format
    YcenterType <- "no"
    YscaleType <- "no"
    if (preProcY != "no") {
        YcenterType <- "mc"
        if (preProcY != "mc")
            YscaleType <- preProcY
    }
    
    # ----- Parameters init
    release <- ""

    YhatDaSave <- list()
    cvTestIndex <- list()
    cvTrainingIndex <- list()
    modelMain <- list()
    
    if (verbose) {
        cat("Please wait... The cross-validation process begins.")
        flush.console()
    }
    
    AllYhat <- c()

    nrun <- if (cvType=='nfold') nfold else nMC
        
    for (icv in 0:nrun) {
        # Update progression bar
        if (verbose) {
            cat("\r", "                                           ", "\r")
            progress <- round(icv * 100 / nrun, 0)
            cat(sprintf("Progression : %.2f%% \r", progress))
            flush.console()
        }
        
        # Set up Cross-Validation
        cvSet <- koplsCrossValSet(K = K, Y = Y, cvFrac = cvFrac, cvType = cvType, 
                                  nfold = nrun, nfoldRound = icv)#, random.seed = 10403+icv)
        #TODO: check when nfold > 1, i.e. repeated test samples
        cvTestIndex[[length(cvTestIndex)+1]] <- cvSet$testIndex
        cvTrainingIndex[[length(cvTrainingIndex)+1]] <- cvSet$trainingIndex
        
        # Get Kernel matrices 
        # % change so that this is done in the K matrix only once and 
        # selected by indices.
        KtrTr <- cvSet$KTrTr
        KteTe <- cvSet$KTeTe
        KteTr <- cvSet$KTeTr
        
        # Center Y and kernel matrices
        YScaleObj <- koplsScale(X = cvSet$yTraining, 
                                centerType = YcenterType,
                                scaleType = YscaleType)
        YScaleObjTest <- koplsScaleApply(model = YScaleObj, X = cvSet$yTest)

        if (preProcK == "mc") {
            KteTe <- koplsCenterKTeTe(KteTe = KteTe, KteTr = KteTr, KtrTr = KtrTr)
            KteTr <- koplsCenterKTeTr(KteTr = KteTr, KtrTr = KtrTr)
            KtrTr <- koplsCenterKTrTr(K = KtrTr)
        }

        # Estimate K-OPLS model
        model <- koplsModel(K = KtrTr, Y = YScaleObj$X, A = maxPcomp, 
                            nox = maxOcomp, preProcK = "no", preProcY = "no")
        
        # Set up model stats
        ssy     <- sum(YScaleObjTest$X^2)
        ssyVars <- sum(YScaleObjTest$X^2)
        #ssx     <- sum(diag(KteTe))
        
        if (icv == 1) {
            ssyTot     <- ssy
            ssyVarsTot <- ssyVars
            #ssxTot     <- ssx
        } else if (icv > 1) {
            ssyTot     <- ssyTot+ssy
            ssyVarsTot <- ssyVarsTot+ssyVars        
            #ssxTot     <- ssxTot+ssx
        }

        modelPredys <- mclapply(1:(maxOcomp+1), mc.cores=mc.cores, function(nOcomp) {
            # Consensus OPLS predict Yhat
            modelPredy <- koplsPredict(KteTr = KteTr, Ktest = KteTe, Ktrain = KtrTr,
                                       model = model, nox = nOcomp-1, 
                                       rescaleY = FALSE)
            return (list(modelPredy=modelPredy,
                         nOcomp=nOcomp))
        })
        AllYhatind <- do.call(cbind, lapply(modelPredys, function(mp) {
            tmp <- koplsRescale(scaleS = YScaleObj, varargin = mp$modelPredy$Yhat)
            colnames(tmp$X) <- paste0(colnames(tmp$X), ifelse(mp$nOcomp==1, "_p", paste0("_po", mp$nOcomp-1)))
            return (tmp$X)
        }))
        
        if (icv==0) 
            pressyTot <-setNames(rep(0,  maxOcomp+1), c("p", paste0("po", 1:maxOcomp))) 
            #pressyTot <- 0 
        else
            pressyTot <- pressyTot + sapply(modelPredys, function(mp) {
                sum((YScaleObjTest$X - mp$modelPredy$Yhat)^2)
            })
        # pressyVarsTot <- sapply(modelPredys, function(mp) {
        #     colSums((YScaleObjTest$X - mp$modelPredy$Yhat)^2) # TODO: clarify: colSums is for Y with variable per column, not for dummy
        # })
        
        # Save Yhat for all rounds
        if (icv == 0) {
            YhatDaSave <- replicate(maxOcomp+1, list())
        }
        
        # + mean on Yhat
        for (mp in modelPredys) { #TODO: reform to lapply
            tmp <- koplsRescale(scaleS = YScaleObj, varargin = mp$modelPredy$Yhat)
            YhatDaSave[[mp$nOcomp]][[length(YhatDaSave[[mp$nOcomp]])+1]] <- tmp$X
            # If highest number of oscs, save Yhat and Xhat
            if (mp$nOcomp == maxOcomp+1) {
                if (icv == 0) {
                    Yhat <- list()
                }
                tmp <- koplsRescale(scaleS = YScaleObj, varargin = mp$modelPredy$Yhat)
                Yhat[[length(Yhat)+1]] <- tmp$X
            }
        }
        
        AllYhat <- rbind(AllYhat, AllYhatind) #TODO: use list
    } # end icv

    if (verbose) {
        cat("Cross-validation is complete.                              ")
        flush.console()
    }
    
    #if (F){
    if (modelType == "da") {
        ## CV DA metrics
        daMetrics_list <- list()
        daMetrics_list$sens <- daMetrics_list$spec <- matrix(nrow=nclasses, ncol=0)
        daMetrics_list$observed <- matrix(classVect[unlist(cvTestIndex[-1])], ncol=1)
        daMetrics_list$predicted <- matrix(NA, 
                                           nrow=nrow(daMetrics_list$observed), 
                                           ncol=maxOcomp+1,
                                           dimnames=list(NULL, c("p", paste0("po", 1:maxOcomp))))
        # Get sens/spec for each y-orth component
        for (i in 1:(maxOcomp + 1)) {
            # Predicted class on test
            predClass <- koplsMaxClassify(X = do.call(rbind, YhatDaSave[[i]][-1]))

            # Calculate sensitivity and specificity
            daMetrics <- koplsSensSpec(trueClass = classVect[unlist(cvTestIndex[-1])],
                                       predClass = predClass,
                                       labelClass = colnames(Y)) 
            daMetrics_list$sens <- cbind(daMetrics_list$sens, 
                                         matrix(daMetrics[, "sens"], 
                                                dimnames=list(rownames(daMetrics), 
                                                              ifelse(i==1, "p", paste0("po", i-1)))))
            daMetrics_list$spec <- cbind(daMetrics_list$spec, 
                                         matrix(daMetrics[, "spec"], 
                                                dimnames=list(rownames(daMetrics), 
                                                              ifelse(i==1, "p", paste0("po", i-1)))))
            #daMetrics_list$tot_sens[i] <- daMetrics[nrow(daMetrics), "sens"]
            #daMetrics_list$meanSens[i] <- daMetrics[nrow(daMetrics), "meanSens"]
            #daMetrics_list$meanSpec[i] <- daMetrics[nrow(daMetrics), "meanSpec"]
            daMetrics_list$confusMatrix[[i]] <- koplsConfusionMatrix(trueClass = classVect[unlist(cvTestIndex[-1])], 
                                                                     predClass = predClass)
            daMetrics_list$predicted[,i] <- predClass
        }
        daMetrics_list$nclasses <- nclasses
        modelMain$cv$da <- daMetrics_list
        
        ## Reprediction DA metrics
        daMetrics_list <- list()
        daMetrics_list$sens <- daMetrics_list$spec <- matrix(nrow=nclasses, ncol=0)
        daMetrics_list$observed <- matrix(classVect[cvTestIndex[[1]]], ncol=1)
        daMetrics_list$predicted <- matrix(NA, nrow=nrow(daMetrics_list$observed), ncol=maxOcomp+1, dimnames=list(NULL, c("p", paste0("po", 1:maxOcomp))))
        # Get sens/spec for each y-orth component
        for (i in 1:(maxOcomp + 1)) {
            # Predicted class on all samples
            predClass <- koplsMaxClassify(X = YhatDaSave[[i]][[1]])
            # Calculate sensitivity and specificity
            daMetrics <- koplsSensSpec(trueClass = classVect[cvTestIndex[[1]]],
                                       predClass = predClass,
                                       labelClass = colnames(Y)) 
            daMetrics_list$sens <- cbind(daMetrics_list$sens, 
                                         matrix(daMetrics[, "sens"], 
                                                dimnames=list(rownames(daMetrics),
                                                              ifelse(i==1, "p", paste0("po", i-1)))))
            daMetrics_list$spec <- cbind(daMetrics_list$spec, 
                                         matrix(daMetrics[, "spec"], 
                                                dimnames=list(rownames(daMetrics),
                                                              ifelse(i==1, "p", paste0("po", i-1)))))
            daMetrics_list$confusMatrix[[i]] <- koplsConfusionMatrix(trueClass = classVect[cvTestIndex[[1]]], 
                                                                     predClass = predClass)
            daMetrics_list$predicted[,i] <- predClass
        }
        daMetrics_list$nclasses <- nclasses
        modelMain$da <- daMetrics_list
        #modelMain$cv$da$args$maxOcomp <- maxOcomp
        #modelMain$cv$da$args$maxPcomp <- maxPcomp
    } else if (modelType == "reg") {
        ## CV metrics
        Metrics_list <- list()
        Metrics_list$observed <- Y[unlist(cvTestIndex[-1]),, drop=F]
        Metrics_list$predicted <- matrix(NA, 
                                         nrow=nrow(Metrics_list$observed), 
                                         ncol=maxOcomp+1, 
                                         dimnames=list(NULL, c("p", paste0("po", 1:maxOcomp))))
        # Get sens/spec for each y-orth component
        for (i in 1:(maxOcomp + 1)) {
            # Predicted values
            Metrics_list$predicted[,i] <- do.call(rbind, YhatDaSave[[i]][-1])
        }
        modelMain$cv$reg <- Metrics_list
        
        ## Reprediction metrics
        Metrics_list <- list()
        Metrics_list$observed <- Y
        Metrics_list$predicted <- matrix(NA,
                                         nrow=nrow(Metrics_list$observed), 
                                         ncol=maxOcomp+1, 
                                         dimnames=list(NULL, c("p", paste0("po", 1:maxOcomp))))
        # Get sens/spec for each y-orth component
        for (i in 1:(maxOcomp + 1)) {
            # Predicted values
            Metrics_list$predicted[,i] <- YhatDaSave[[i]][[1]]
        }
        modelMain$reg <- Metrics_list
    }
    #}
    modelMain$Model <- koplsModel(K = K, Y = Y, A = maxPcomp, nox = maxOcomp, 
                                  preProcK = preProcK, preProcY = preProcY)

    #modelMain$cv$Yhat <- do.call(rbind, Yhat[-1])
    #colnames(modelMain$cv$Yhat) <- colnames(Y)
    modelMain$cv$AllYhat <- AllYhat[-(1:nrow(Y)),,drop=F]
    #modelMain$cv$Tcv     <- modelMain$cv$Yhat %*% modelMain$Model$Cp %*% modelMain$Model$Bt[[maxOcomp + 1]] #TODO: why Bt[[maxOcomp + 1]]?
    modelMain$cv$Q2Yhat     <- 1 - pressyTot/ssyTot
    #modelMain$cv$Q2YhatVars <- 1 - pressyVarsTot/ssyVarsTot
    
    modelMain$cv$cvTestIndex     <- cvTestIndex[-1]
    #modelMain$cv$cvTrainingIndex <- cvTrainingIndex[-1]
    
    #modelMain$class <- "koplscv"
    
    return (modelMain)
}


