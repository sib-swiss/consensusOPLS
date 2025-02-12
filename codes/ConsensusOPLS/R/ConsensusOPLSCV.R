#' @title ConsensusOPLSCV
#' @description Function for performing Kernel-OPLS cross-validation for a set 
#' of Y-orthogonal components.
#' This function was first implemented in KOPLS1.1 MATLAB package, and was
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
#' for mean-centering, \code{uv} for mean-centering and scaling to
#' unit-variance, \code{pa} for mean-centering and Pareto-scaling or \code{no}
#' for no mean-centering and no scaling.
#' @param cvFrac numeric. Fraction of observations in the training set during
#' cross-validation. 
#' @param modelType character. Type of model used for the ConsensusOPLS method, 
#' \code{da} for discriminant analysis or \code{reg} for regression. 
#' If \code{da} (default), sensitivity and specificity will be calculated.
#' @param verbose logical. Indicates whether the user wants to see the 
#' progress bar printed. If \code{FALSE}, no output will be printed. If 
#' \code{TRUE} (default) some output will be printed regarding the 
#' cross-validation progress.
#'
#' @returns A list with the following entries (diagnostic parameters which can
#' be used to determine the optimal number of model components):
#' \item{cv}{ list. The cross-validation results.}
#'      \item{Yhat}{ matrix. Predicted Y values.}
#'      \item{AllYhat}{ matrix. All predicted Y values.}
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
#'      \item{nclass}{ integer. Number of classes in model.}
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
    if (!is.matrix(K))
        stop("K is not a matrix.")
    if (!is.numeric(maxPcomp))
        stop("maxPcomp is not numeric.")
    if (!is.numeric(maxOcomp))
        stop("maxOcomp is not numeric.")
    if (!is.numeric(nfold))
        stop("nfold is not numeric.")
    if (!is.numeric(nMC))
        stop("nfold is not numeric.")
    if (!is.numeric(cvFrac))
        stop("cvFrac is not numeric.")
    if (!is.logical(verbose))
        stop("verbose must be a logical value.")
    
    cvType <- match.arg(cvType,
                        choices = c("nfold", "mccv", "mccvb"),
                        several.ok = F)
    preProcK <- match.arg(preProcK,
                          choices = c("mc", "no"),
                          several.ok = F)
    preProcY <- match.arg(preProcY,
                          choices = c("mc", "no", "uv", "pa"),
                          several.ok = F)
    modelType <- match.arg(modelType,
                           choices = c("da", "reg"),
                           several.ok = F)
    
    # ----- Default values control
    if (cvType == "mccvb" && modelType != "da")
        stop("Class balanced MC CV only applicable to `da` modelling.")

    # ----- Variable format control (part 2)
    if (modelType == "da") {
        # Check the response matrix
        if (! (all(Y %in% c(0, 1)) && ncol(Y) > 1))
            stop("Y is not a dummy matrix.")
        
        classVect <- koplsReDummy(Y = Y)
        nclass <- length(unique(x = classVect))
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
        cvSet <- koplsCrossValSet(K = K, Y = Y,
                                  cvFrac = cvFrac, cvType = cvType,
                                  nfold = nrun, nfoldRound = icv)
        #, random.seed = 10403+icv)
        # TODO: check when nfold > 1, i.e. repeated test samples
        cvTestIndex[[length(cvTestIndex)+1]] <- cvSet$testIndex
        cvTrainingIndex[[length(cvTrainingIndex)+1]] <- cvSet$trainingIndex
        
        # Get Kernel matrices 
        # TODO: change so that this is done in the K matrix only once and 
        # selected by indices.
        KtrTr <- cvSet$KTrTr
        KteTe <- cvSet$KTeTe
        KteTr <- cvSet$KTeTr
        
        # Scale Y
        YScaleObj <- koplsScale(X = cvSet$yTraining,
                                centerType = YcenterType,
                                scaleType = YscaleType)
        YScaleObjTest <- koplsScaleApply(model = YScaleObj,
                                         X = cvSet$yTest)

        # Center kernel matrices
        if (preProcK == "mc") {
            KteTe <- koplsCenterKTeTe(KteTe = KteTe,
                                      KteTr = KteTr,
                                      KtrTr = KtrTr)
            KteTr <- koplsCenterKTeTr(KteTr = KteTr,
                                      KtrTr = KtrTr)
            KtrTr <- koplsCenterKTrTr(KtrTr = KtrTr)
        }

        # Estimate K-OPLS model
        model <- koplsModel(K = KtrTr, Y = YScaleObj$X,
                            A = maxPcomp, nox = maxOcomp,
                            preProcK = "no", preProcY = "no")
        
        # Set up model stats
        ssy     <- sum(YScaleObjTest$X^2)
        #ssyVars <- sum(YScaleObjTest$X^2)
        #ssx     <- sum(diag(KteTe))
        
        if (icv == 1) {
            ssyTot     <- ssy
            #ssyVarsTot <- ssyVars
            #ssxTot     <- ssx
        } else if (icv > 1) {
            ssyTot     <- ssyTot+ssy
            #ssyVarsTot <- ssyVarsTot+ssyVars        
            #ssxTot     <- ssxTot+ssx
        }

        modelPredys <- mclapply(
            1:(maxOcomp+1), mc.cores=mc.cores, function(nOcomp) {
                # predict Yhat with (nOcomp-1) orthogonal components
                modelPredy <- koplsPredict(KteTr = KteTr,
                                           Ktest = KteTe,
                                           Ktrain = KtrTr,
                                           model = model,
                                           nox = nOcomp-1,
                                           rescaleY = FALSE)
                return (list(modelPredy = modelPredy,
                             nOcomp = nOcomp))
            })
        
        AllYhatind <- do.call(cbind, lapply(modelPredys, function(mp) {
            tmp <- koplsRescale(scaleS = YScaleObj,
                                varargin = mp$modelPredy$Yhat)
            colnames(tmp$X) <- paste0(colnames(tmp$X),
                                      ifelse(mp$nOcomp==1,
                                             "_p",
                                             paste0("_po", mp$nOcomp-1)))
            return (tmp$X)
        }))
        
        if (icv==0) {
            pressyTot <- setNames(rep(0,  maxOcomp+1),
                                  c("p", paste0("po", 1:maxOcomp)))
        } else {
            pressyTot <- pressyTot + sapply(modelPredys, function(mp) {
                sum((YScaleObjTest$X - mp$modelPredy$Yhat)^2)
            })
        }
        # pressyVarsTot <- sapply(modelPredys, function(mp) {
        #     colSums((YScaleObjTest$X - mp$modelPredy$Yhat)^2)
        # TODO: clarify: colSums is for Y with variable per column, not for dummy
        # })
        
        # Save Yhat for all rounds
        if (icv == 0) {
            YhatDaSave <- replicate(maxOcomp+1, list())
        }
        
        # Add mean on Yhat by koplsRescale
        for (mp in modelPredys) { #TODO: reform to lapply
            tmp <- koplsRescale(scaleS = YScaleObj,
                                varargin = mp$modelPredy$Yhat)
            YhatDaSave[[mp$nOcomp]][[length(YhatDaSave[[mp$nOcomp]])+1]] <-
                tmp$X
            # # Store Yhat at the highest number of orthogonal components
            # if (mp$nOcomp == maxOcomp+1) {
            #     if (icv == 0) Yhat <- list()
            #     tmp <- koplsRescale(scaleS = YScaleObj,
            #                         varargin = mp$modelPredy$Yhat)
            #     Yhat[[length(Yhat)+1]] <- tmp$X
            # }
        }
        
        AllYhat <- rbind(AllYhat, AllYhatind) #TODO: use list
    } # end icv

    if (verbose) {
        cat("Cross-validation is complete.")
        flush.console()
    }
    
    #if (F){
    if (modelType == "da") {
        ## CV DA metrics
        daMetrics_list <- list()
        daMetrics_list$sens <- daMetrics_list$spec <-
            matrix(nrow=nclass, ncol=0)
        daMetrics_list$observed <- matrix(classVect[unlist(cvTestIndex[-1])],
                                          ncol=1)
        daMetrics_list$predicted <- matrix(
            NA, 
            nrow=nrow(daMetrics_list$observed),
            ncol=maxOcomp+1,
            dimnames=list(NULL, c("p", paste0("po", 1:maxOcomp))))
        # Get sens/spec for each Y-orth component
        for (i in 1:(maxOcomp + 1)) {
            # Predicted class on test
            predClass <- koplsMaxClassify(X = do.call(rbind,
                                                      YhatDaSave[[i]][-1]))

            # Calculate sensitivity and specificity
            daMetrics <- koplsSensSpec(
                trueClass = classVect[unlist(cvTestIndex[-1])],
                predClass = predClass,
                labelClass = colnames(Y)) 
            daMetrics_list$sens <- cbind(
                daMetrics_list$sens, 
                matrix(daMetrics[, "sens"], 
                       dimnames=list(rownames(daMetrics), 
                                     ifelse(i==1, "p", paste0("po", i-1)))))
            daMetrics_list$spec <- cbind(
                daMetrics_list$spec, 
                matrix(daMetrics[, "spec"], 
                       dimnames=list(rownames(daMetrics), 
                                     ifelse(i==1, "p", paste0("po", i-1)))))
            #daMetrics_list$tot_sens[i] <- daMetrics[nrow(daMetrics), "sens"]
            #daMetrics_list$meanSens[i] <- daMetrics[nrow(daMetrics), "meanSens"]
            #daMetrics_list$meanSpec[i] <- daMetrics[nrow(daMetrics), "meanSpec"]
            daMetrics_list$confusMatrix[[i]] <- koplsConfusionMatrix(
                trueClass = classVect[unlist(cvTestIndex[-1])], 
                predClass = predClass)
            daMetrics_list$predicted[,i] <- predClass
        }
        daMetrics_list$nclass <- nclass
        modelMain$cv$da <- daMetrics_list
        
        ## Reprediction DA metrics
        daMetrics_list <- list()
        daMetrics_list$sens <- daMetrics_list$spec <-
            matrix(nrow=nclass, ncol=0)
        daMetrics_list$observed <- matrix(classVect[cvTestIndex[[1]]], ncol=1)
        daMetrics_list$predicted <- matrix(
            NA,
            nrow=nrow(daMetrics_list$observed),
            ncol=maxOcomp+1,
            dimnames=list(NULL, c("p", paste0("po", 1:maxOcomp))))
        # Get sens/spec for each y-orth component
        for (i in 1:(maxOcomp + 1)) {
            # Predicted class on all samples
            predClass <- koplsMaxClassify(X = YhatDaSave[[i]][[1]])
            # Calculate sensitivity and specificity
            daMetrics <- koplsSensSpec(trueClass = classVect[cvTestIndex[[1]]],
                                       predClass = predClass,
                                       labelClass = colnames(Y)) 
            daMetrics_list$sens <- cbind(
                daMetrics_list$sens,
                matrix(daMetrics[, "sens"], 
                       dimnames=list(rownames(daMetrics),
                                     ifelse(i==1, "p", paste0("po", i-1)))))
            daMetrics_list$spec <- cbind(
                daMetrics_list$spec, 
                matrix(daMetrics[, "spec"], 
                       dimnames=list(rownames(daMetrics),
                                     ifelse(i==1, "p", paste0("po", i-1)))))
            daMetrics_list$confusMatrix[[i]] <- koplsConfusionMatrix(
                trueClass = classVect[cvTestIndex[[1]]], 
                predClass = predClass)
            daMetrics_list$predicted[,i] <- predClass
        }
        daMetrics_list$nclass <- nclass
        modelMain$cv$da <- daMetrics_list
    } else if (modelType == "reg") {
        ## CV regression metrics
        Metrics_list <- list()
        Metrics_list$observed <- Y[unlist(cvTestIndex[-1]), , drop=F]
        Metrics_list$predicted <- matrix(
            NA, 
            nrow=nrow(Metrics_list$observed), 
            ncol=maxOcomp+1, 
            dimnames=list(NULL, c("p", paste0("po", 1:maxOcomp))))
        # Get sens/spec for each Y-orth component
        for (i in 1:(maxOcomp + 1)) {
            # Predicted values
            Metrics_list$predicted[, i] <- do.call(rbind, YhatDaSave[[i]][-1])
        }
        modelMain$cv$reg <- Metrics_list
        
        ## Reprediction metrics
        Metrics_list <- list()
        Metrics_list$observed <- Y
        Metrics_list$predicted <- matrix(
            NA,
            nrow=nrow(Metrics_list$observed), 
            ncol=maxOcomp+1, 
            dimnames=list(NULL, c("p", paste0("po", 1:maxOcomp))))
        # Get sens/spec for each Y-orth component
        for (i in 1:(maxOcomp + 1)) {
            # Predicted values
            Metrics_list$predicted[,i] <- YhatDaSave[[i]][[1]]
        }
        modelMain$reg <- Metrics_list
    }
    #}
    
    # modelMain$Model <- koplsModel(K = K, Y = Y, A = maxPcomp, nox = maxOcomp, 
    #                               preProcK = preProcK, preProcY = preProcY)

    #modelMain$cv$Yhat <- do.call(rbind, Yhat[-1])
    #colnames(modelMain$cv$Yhat) <- colnames(Y)
    modelMain$cv$AllYhat <- AllYhat[-(1:nrow(Y)), , drop=F]
    #modelMain$cv$Tcv     <- modelMain$cv$Yhat %*% modelMain$Model$Cp %*% modelMain$Model$Bt[[maxOcomp + 1]] #TODO: why Bt[[maxOcomp + 1]]?
    modelMain$cv$Q2Yhat <- 1 - pressyTot/ssyTot
    #modelMain$cv$Q2YhatVars <- 1 - pressyVarsTot/ssyVarsTot
    
    modelMain$cv$cvTestIndex <- cvTestIndex[-1]
    #modelMain$cv$cvTrainingIndex <- cvTrainingIndex[-1]
    
    class(modelMain) <- "koplscv"
    
    return (modelMain)
}


