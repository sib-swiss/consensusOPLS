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
#' @param A numeric. The number of Y-predictive components (integer). 
#' @param oax numeric. The number of Y-orthogonal components (integer).
#' @param nbrcv numeric. Number of cross-validation rounds (integer).
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
#' cross-validation (integer). 
#' @param modelType: character. Type of model used for the ConsensusOPLS method. 
#' It could be \code{da} for discriminant analysis or \code{reg} for regression. 
#' If \code{da} (default), sensitivity and specificity will be calculated.
#' @param verbose logical. Indicates whether the user wants to see the 
#' progress bar printed. If \code{FALSE}, no output will be printed. If 
#' \code{TRUE} (default) some output will be printed regarding the 
#' cross-validation progress.
#'
#' @return A list with the following entries (diagnostic parameters which can be 
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
#'            \item{oax}{ integer. Number of Y-orthogonal components.}
#'            \item{A}{ integer. Number of Y-predictive components.}
#' \item{class}{ character. Model class is \code{koplscv}.}
#' @examples
#' #TO DO
#' @importFrom utils flush.console
#' @keywords internal  

ConsensusOPLSCV <- function(K, Y, 
                            A, oax, nbrcv, 
                            cvType,
                            preProcK = "no", preProcY = "no", 
                            cvFrac, 
                            modelType = "da", verbose = TRUE){
    # ----- Variable format control (part 1)
    if (!is.matrix(K)) stop("K is not a matrix.")
    if (!is.matrix(Y)) stop("Y is not a matrix.")
    if (!is.numeric(A)) stop("A is not numeric.")
    if (!is.numeric(oax)) stop("oax is not numeric.")
    if (!is.numeric(nbrcv)) stop("nbrcv is not numeric.")
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
    else if(!(modelType %in% c("da", "re"))) stop("modelType must me `da` or `re`.")
    if (!is.logical(verbose)) stop("verbose must be `TRUE` or `FALSE`.")
    
    # ----- Default values control
    if (cvType == "mccvb" & modelType != "da")
        stop("Class balanced MC CV only applicable to `da` modelling.")
    print("conoplscv")
    print(modelType)
    print(head(Y))
    # ----- Variable format control (part 2)
    if (modelType == "da") {
        # Define a parameter for DA decision rule
        drRule <- "max"
        
        # Check the response matrix
      #   if (all(Y %in% c(0, 1))) {
      #       if (ncol(Y) == 1) 
      #           Y <- koplsDummy(X = Y, numClasses = NA)
      #       classVect <- koplsReDummy(Y = Y)
      #   } else {
      #       if (ncol(Y) == 1) {
      #           classVect <- Y
      #           Y <- koplsDummy(X = Y+1, numClasses = NA)
      #       } else
      #           stop("modelType is `da`, but Y appears to be neither dummy matrix nor 
      # a vector of (integer) class labels.")
      #   }
        classVect <- as.vector(Y)
        nclasses <- length(unique(x = classVect))
    }

    # ----- Convert Y-scaling to more explicit format
    YcenterType <- "no"
    YscaleType <- "no"
    if (preProcY != "no") {
        YcenterType <- "mc"
        if(preProcY != "mc")
            YscaleType <- preProcY
        
    }
    
    # ----- Parameters init
    release <- ""

    pressy <- matrix(data = 0, nrow = oax+1, ncol = 1)
    pressyVars <- matrix(data = 0, nrow = oax+1, ncol = 1)
    pressyTot <- matrix(data = 0, nrow = oax+1, ncol = 1)
    pressyVarsTot <- matrix(data = 0, nrow = oax+1, ncol = 1)
    YhatDaSave <- list()
    cvTestIndex <- c()
    cvTrainingIndex <- c()
    
    if (verbose) {
        cat("Please wait... The cross-validation process begins.")
        flush.console()
    }
    
    AllYhat <- c()
    
    for (icv in 1:nbrcv) {
        # Update progression bar
        if (verbose) {
            cat("\r", "                                           ", "\r")
            progress <- round(icv * 100 / nbrcv, 0)
            cat(sprintf("Progression : %.2f%% \r", progress))
            flush.console()
        }
        
        # Set up Cross-Validation
        cvSet <- koplsCrossValSet(K = K, Y = Y, modelFrac = cvFrac, type = cvType, 
                                  nfold = nbrcv, nfoldRound = icv)
        cvTestIndex <- c(cvTestIndex, cvSet$testIndex)
        cvTrainingIndex <- c(cvTrainingIndex, cvSet$trainingIndex)
        
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
        model <- koplsModel(K = KtrTr, Y = YScaleObj$matrix, A = A, 
                            nox = oax, preProcK = "no", preProcY = "no")
        
        # Set up model stats
        ssy <- sum(YScaleObjTest$matrix^2)
        ssyVars <- sum(YScaleObjTest$matrix^2)
        ssx <- sum(diag(KteTe))
        
        if (icv == 1) {
            ssyTot <- ssy
            ssyVarsTot <- ssyVars
            ssxTot <- ssx
        } else {
            ssyTot <- ssyTot+ssy
            ssyVarsTot <- ssyVarsTot+ssyVars        
            ssxTot <- ssxTot+ssx
        }
        
        # for each combination of Y-osc components
        AllYhatind <- c()
        
        for (ioax in 1:(oax+1)) {
            ioay <- 1
            # Consensus OPLS predict Yhat
            modelPredy <- koplsPredict(KteTr = KteTr, Ktest = KteTe, Ktrain = KtrTr,
                                       model = model, nox = ioax-1, 
                                       rescaleY = FALSE)
            tmp <- koplsRescale(scaleS = YScaleObj, varargin = modelPredy$Yhat)
            AllYhatind <- cbind(AllYhatind, tmp$X)
            pressy[ioax, ioay] <- sum((YScaleObjTest$matrix - 
                                           modelPredy$Yhat)^2)
            pressyVars[ioax, ioay] <- sum((YScaleObjTest$matrix - 
                                               modelPredy$Yhat)^2)
            
            if (icv == 1) {
                pressyTot[ioax, ioay] <- pressy[ioax, ioay]
                pressyVarsTot[ioax, ioay] <- pressyVars[ioax, ioay]
            } else {
                pressyTot[ioax,ioay] <- pressyTot[ioax,ioay] + pressy[ioax,ioay]
                pressyVarsTot[ioax,ioay] <- pressyVarsTot[ioax,ioay] + pressyVars[ioax,ioay]
            }
            
            # If 'da', save Yhat for all rounds
            if (modelType == "da") {
                if (icv == 1) {
                    YhatDaSave[[ioax]] <- list()
                }
                
                # + mean on Yhat
                tmp <- koplsRescale(YScaleObj, modelPredy$Yhat)
                YhatDaSave[[ioax]] <- rbind(YhatDaSave[[ioax]], tmp$X)
            }
            
            # If highest number of oscs, save Yhat and Xhat
            if (ioax == oax+1) {  # && ioay == oay + 1) 
                if (icv == 1) {
                    Yhat <- list()
                }
                tmp <- koplsRescale(YScaleObj, modelPredy$Yhat)
                Yhat <- rbind(Yhat, tmp$X)
            }
        }
        AllYhat <- rbind(AllYhat, AllYhatind)
    } # end icv
    
    if (verbose) {
        cat("Cross-validation is complete.                              ")
        flush.console()
    }
    
    KtrTr <- K
    modelMain <- list()
    modelMain$Model <- koplsModel(K = KtrTr, Y = Y, A = A, nox = oax, 
                                  preProcK = preProcK, preProcY = preProcY)
    modelMain$cv$Yhat <- matrix(data = unlist(Yhat), 
                                ncol = length(Yhat[1,]))
    modelMain$cv$AllYhat <- AllYhat
    modelMain$cv$Tcv <- modelMain$cv$Yhat %*% modelMain$Model$Cp %*% modelMain$Model$Bt[[oax + 1]]
    modelMain$cv$Q2Yhat <- 1 - pressyTot/ssyTot
    modelMain$cv$Q2YhatVars <- 1 - pressyVarsTot/ssyVarsTot
    modelMain$cv$cvTestIndex <- cvTestIndex
    modelMain$cv$cvTrainingIndex <- cvTrainingIndex
    
    print("debug")
    print(head(classVect))
    daMetrics_list <- list()
    if (modelType == "da") {
        # Get sens/spec for each y-orth component
        for (i in 1:(oax + 1)) {
            if (drRule == "max") {
                predClass <- koplsMaxClassify(X = YhatDaSave[[i]])
            } else if (drRule == "fixed") {
                predClass <- koplsBasicClassify(X = YhatDaSave[[i]], 
                                                k = 1/nclasses)
            } else {
                warning(paste0('Decision rule given: ', drRule, 
                               ' is not valid/implemented.'))
            }
            
            # Calculate sensitivity and specificity
            daMetrics <- koplsSensSpec(trueClass = classVect[cvTestIndex],
                                       predClass = predClass)
            daMetrics_list$sens[i] <- daMetrics[i, "sens"]
            daMetrics_list$spec[i] <- daMetrics[i, "spec"]
            daMetrics_list$tot_sens[i] <- daMetrics[nrow(daMetrics), "sens"]
            daMetrics_list$meanSens[i] <- daMetrics[nrow(daMetrics), "meanSens"]
            daMetrics_list$meanSpec[i] <- daMetrics[nrow(daMetrics), "meanSpec"]
        }
        
        # Calculate sensitivity and specificity
        daMetrics_list$confusMatrix <- koplsConfusionMatrix(trueClass = classVect[cvTestIndex], 
                                                            predClass = predClass)
        daMetrics_list$trueClass <- classVect[cvTestIndex]
        daMetrics_list$nclasses <- nclasses
        modelMain$da <- daMetrics_list
        modelMain$da$predClass <- predClass
        modelMain$da$decisionRule <- drRule
        
        # Change to original order if NFOLD CV
        if (cvType == "nfold") {
            cvOrder <- sort(x = cvTestIndex, decreasing = FALSE)
            modelMain$da$predClass <- modelMain$da$predClass[cvOrder]
            modelMain$da$trueClass <- modelMain$da$trueClass[cvOrder]
        }
    }
    
    # Change to original order if NFOLD CV
    if (cvType == "nfold") {
        cvOrder <- sort(x = cvTestIndex, decreasing = FALSE)
        modelMain$cv$Yhat <- modelMain$cv$Yhat[cvOrder, ]
        modelMain$cv$Tcv <- modelMain$cv$Tcv[cvOrder, ]
    }
    
    modelMain$class <- "koplscv"
    modelMain$da$args$oax <- oax
    modelMain$da$args$A <- A
    
    return (modelMain)
}


