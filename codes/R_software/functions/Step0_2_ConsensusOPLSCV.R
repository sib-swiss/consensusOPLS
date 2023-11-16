#' ConsensusOPLSCV
#' Function for performing Kernel-OPLS cross-validation for a set of 
#' Y-orthogonal components.
#' This function was first implemented in KOPLS1.1 MATLAB package, but was
#' modified in 2012 for the ConsensusOPLS method.
#'
#' @param K: The kernel matrix (un-centered); see 'koplsKernel()' for details.
#' @param Y: The response matrix (un-centered/scaled). Could be binary (for
#' discriminant analysis) or real-valued (for classical OPLS analysis).
#' @param A: The number of Y-predictive components (integer). 
#' @param oax: The number of Y-orthogonal components (integer).
#' @param nbrcv: Number of cross-validation rounds (integer).
#' @param cvType: Type of cross-validation used. Either `nfold` for n-fold
#' cross-validation, `mccv` for Monte Carlo CV or `mccvb` for Monte Carlo 
#' class-balanced CV. See also 'koplsCrossValSet()' for details. 
#' @param preProcK: Pre-processing settings for the kernel matrix. Either 
#' `mc` for mean-centering or `no` for no pre-processing.
#' @param preProcY: Pre-processing parameter for Y. Either `mc` for 
#' mean-centering, `uv` for mc + scaling to unit-variance, `pareto` for 
#' mc + Pareto-scaling or `no` for no scaling.
#' @param cvFrac: Fraction of observations in the training set during 
#' cross-validation (integer). 
#' @param modelType: character which indicates the type of model used for the 
#' ConsensusOPLS method. It could be `da` for discriminant analysis or `re` 
#' for regression. If 'da' (default), sensitivity and specificity will be 
#' calculated.
#' @param verbose: logical which indicates whether the user wants to see the 
#' progress bar printlayed. If `FALSE`, no output will be printlayed. If `TRUE` 
#' (default) some output will be printlayed regarding the cross-validation 
#' progress .
#'
#' @return model: a number of diagnostic parameters which can be used to 
#' determine the optimal number of model components. It contains :
#' `Model`: list. The training a K-OPLS model.
#'      `Cp`: matrix. Y loading matrix.
#'      `Sp`: matrix. Sigma matrix, containing singular values from Y'*K*Y used 
#'      for scaling.
#'      `Sps`: matrix. Scaled Sigma matrix, containing scaled singular values.
#'      `Up`: matrix. Y score matrix.
#'      `Tp`: list. Predictive score matrix for all Y-orthogonal components.
#'      `T`: matrix. Predictive score matrix for the final model.
#'      `co`: list. Y-orthogonal loading vectors.
#'      `so`: list. Eigenvalues from estimation of Y-orthogonal loading vectors.
#'      `to`: list. Weight vector for the i-th latent component of the KOPLS 
#'      model.
#'      `To`: matrix. Y-orthogonal score matrix.
#'      `toNorm`: list. Norm of the Y-orthogonal score matrix prior to scaling.
#'      `Bt`: list. T-U regression coefficients for predictions.
#'      `A`: numeric. Number of predictive components.
#'      `nox`: numeric. Number of Y-orthogonal components.
#'      `K`: matrix. The kernel matrix.
#'      `EEprime`: matrix. The deflated kernel matrix for residual statistics.
#'      `sstot_K`: numeric. Total sums of squares in 'K'.
#'      `R2X`: numeric. Cumulative explained variation for all model components.
#'      `R2XO`: numeric. Cumulative explained variation for Y-orthogonal model 
#'      components.
#'      `R2XC`: numeric. Explained variation for predictive model components 
#'      after addition of Y-orthogonal model components.
#'      `sstot_Y`: numeric. Total sums of squares in Y.
#'      `R2Y`: numeric. Explained variation of Y.
#'      `R2Yhat`: numeric. Variance explained by the i-th latent component of 
#'      the model.
#'      `preProc$K`: character. Pre-processing setting for K.
#'      `preProc$Y`: character. Pre-processing setting for Y.
#'      `preProc$paramsY`: character. Pre-processing scaling parameters for Y.
#' `cv`: list. The cross-validation results.
#'      `Yhat`: matrix. Predicted Y values .
#'      `AllYhat`: matrix. All predicted Y values as a concatenated matrix.
#'      `Tcv`: matrix. redictive score vector T for all cross-validation rounds.
#'      `Q2Yhat`: matrix. Total Q-square result for all Y-orthogonal components.
#'      `Q2YhatVars`: matrix. Q-square result per Y-variable for all 
#'      Y-orthogonal components.
#'      `cvTestIndex`: matrix. Indices for the test set observations during the
#'      cross-validation rounds.
#'      `cvTrainIndex`: matrix. Indices for the training set observations during 
#'      the cross-validation rounds.
#' `da`: list. Cross-validation results specifically for discriminant analysis:
#'      `totalResults `: list. The results over all classes and averages.
#'      `classResults `: list. The results for each class.
#'      `predClass`: integer. Predicted class list per class and Y-orthogonal
#'      components.
#'      `trueClass`: integer. Predicted class list per class and Y-orthogonal 
#'      components.
#'      `sensSpec`: integer. Sensitivity and specificity values per class and
#'      Y-orthogonal components.
#'      `confusionMatrix`: matrix. Confusion matrix during cross-validation
#'      rounds.
#'      `nclasses`: integer. Number of classes in model.
#'      `decisionRule`: character. Decision rule used: 'max' or 'fixed'.
#'      `args`: list. Arguments to the function:
#'            `oax`: integer. Number of Y-orthogonal components.
#'            `A`: integer. Number of Y-predictive components.
#' `class`: character. Model class: koplscv.
#' 
#' @examples
#' #TO DO

ConsensusOPLSCV <- function(K, Y, 
                            A, oax, nbrcv, 
                            cvType,
                            preProcK = "no", preProcY = "no", 
                            cvFrac, 
                            modelType = "da", verbose = TRUE){
  # ----- Variable format control (part 1)
  if(!is.matrix(K)){stop("K is not a matrix.")}
  if(!is.matrix(Y)){stop("Y is not a matrix.")}
  if(!is.numeric(A)){stop("A is not numeric.")}
  if(!is.numeric(oax)){stop("oax is not numeric.")}
  if(!is.numeric(nbrcv)){stop("nbrcv is not numeric.")}
  if(!is.character(cvType)){stop("cvType is not a character.")
  } else{
    if(!(cvType %in% c("nfold", "mccv", "mccvb"))){
      stop("cvType must be `nfold`, `mccv` or `mccvb`.")}
  }
  if(!is.character(preProcK)){stop("preProcK is not a character.")
  } else{
    if(!(preProcK %in% c("mc", "no"))){stop("preProcK must be `mc` or `no`.")}
  }
  if(!is.character(preProcY)){stop("preProcY is not a character.")
  } else{ 
    if(!(preProcY %in% c("mc", "uv", "pareto", "no"))){
      stop("preProcY must be `mc`, `uv`, `pareto` or `no`.")}
  }
  if(!is.numeric(cvFrac)){stop("cvFrac is not numeric.")}
  if(!is.character(modelType)){stop("modelType is not a character.")
  } else{
    if(!(modelType %in% c("da", "re"))){stop("modelType must me `da` or `re`.")}
  }
  if(!is.logical(verbose)){stop("verbose must be `TRUE` or `FALSE`.")}
  
  # ----- Function loading control
  if(!exists("koplsDummy", mode = "function")){
    warning("Remember to load the source code for the `koplsDummy` function.")
  }
  if(!exists("koplsReDummy", mode = "function")){
    warning("Remember to load the source code for the `koplsReDummy` function.")
  }
  if(!exists("koplsCrossValSet", mode = "function")){
    warning("Remember to load the source code for the `koplsCrossValSet` function.")
  }
  if(!exists("koplsScale", mode = "function")){
    warning("Remember to load the source code for the `koplsScale` function.")
  }
  if(!exists("koplsScaleApply", mode = "function")){
    warning("Remember to load the source code for the `koplsScaleApply` function.")
  }
  if(!exists("koplsCenterKTeTe", mode = "function")){
    warning("Remember to load the source code for the `koplsCenterKTeTe` function.")
  }
  if(!exists("koplsCenterKTeTr", mode = "function")){
    warning("Remember to load the source code for the `koplsCenterKTeTr` function.")
  }
  if(!exists("koplsCenterKTrTr", mode = "function")){
    warning("Remember to load the source code for the `koplsCenterKTrTr` function.")
  }
  if(!exists("koplsModel", mode = "function")){
    warning("Remember to load the source code for the `koplsModel` function.")
  }
  if(!exists("koplsPredict", mode = "function")){
    warning("Remember to load the source code for the `koplsPredict` function.")
  }
  if(!exists("koplsRescale", mode = "function")){
    warning("Remember to load the source code for the `koplsRescale` function.")
  }
  if(!exists("koplsMaxClassify", mode = "function")){
    warning("Remember to load the source code for the `koplsMaxClassify` function.")
  }
  if(!exists("koplsBasicClassify", mode = "function")){
    warning("Remember to load the source code for the `koplsBasicClassify` function.")
  }
  if(!exists("koplsSensSpec", mode = "function")){
    warning("Remember to load the source code for the `koplsSensSpec` function.")
  }
  if(!exists("koplsConfusionMatrix", mode = "function")){
    warning("Remember to load the source code for the `koplsConfusionMatrix` function.")
  }

  # ----- Default values control
  if(cvType == "mccvb" & modelType != "da"){
    stop("Class balanced MC CV only applicable to `da` modelling.")
  } 

  # ----- Variable format control (part 2)
  if(modelType == "da"){
    # Define a parameter for DA decision rule
    drRule <- "max"
    
    # Check the response matrix
    temp <- unique(Y)
    if(all(temp == 0 | temp == 1)){
      if(ncol(Y) == 1){Y <- koplsDummy(Y)}
      classVect <- koplsReDummy(Y)
    } else{
      if(all(mod(Y,1)) == 0 & ncol(Y) == 1){
        classVect <- Y
        Y <- koplsDummy(Y+1)
      } else{
      stop("modelType is `da`, but Y appears to be neither dummy matrix nor 
      a vector of (integer) class labels.")}
    }
    nclasses <- base::length(base::unique(classVect))
  }
  
  # ----- Convert Y-scaling to more explicit format
  YcenterType <- "no"
  YscaleType <- "no"
  if(preProcY != "no"){
    YcenterType <- "mc"
    if(preProcY != "mc"){
      YscaleType <- preProcY
    }
  }
  
  # ----- Parameters init
  release <- ""
  set.seed(1214)
  
  pressy <- base::matrix(data = 0, nrow = oax+1, ncol = 1)
  pressyVars <- base::matrix(data = 0, nrow = oax+1, ncol = 1)
  pressyTot <- base::matrix(data = 0, nrow = oax+1, ncol = 1)
  pressyVarsTot <- base::matrix(data = 0, nrow = oax+1, ncol = 1)
  YhatDaSave <- list()
  cvTestIndex <- c()
  cvTrainingIndex <- c()

  if(verbose){
    base::cat("Please wait... The cross-validation process begins.")
    utils::flush.console()
  }
  
  AllYhat <- c()
  
  for(icv in 1:nbrcv){
    # Update progression bar
    if(verbose){
      cat("\r", "                                           ", "\r")
      progress <- round(icv * 100 / nbrcv, 0)
      base::cat(base::sprintf("Progression : %.2f%% \r", progress))
      utils::flush.console()
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
    
    if(preProcK == "mc"){
      KteTe <- koplsCenterKTeTe(KteTe,KteTr,KtrTr)
      KteTr <- koplsCenterKTeTr(KteTr,KtrTr)
      KtrTr <- koplsCenterKTrTr(KtrTr)
    }
    
    # Estimate K-OPLS model
    model <- koplsModel(K = KtrTr, Y = YScaleObj$matrix, A = A, 
                        nox = oax, preProcK = "no", preProcY = "no")
    
    # Set up model stats
    ssy <- base::sum( base::sum(YScaleObjTest$matrix**2 ))
    ssyVars <- base::sum(YScaleObjTest$matrix**2)
    ssx <- base::sum( base::diag(KteTe))
    
    if(icv == 1){
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
    
    for(ioax in 1:(oax+1)){
      for(ioay in 1:1){
        # Consensus OPLS predict Yhat
        modelPredy <- koplsPredict(KteTr = KteTr, Ktest = KteTe, Ktrain = KtrTr,
                                   model = model, nox = ioax-1, 
                                   rescaleY = FALSE)
        tmp <- koplsRescale(scaleS = YScaleObj, varargin = modelPredy$Yhat)
        AllYhatind <- cbind(AllYhatind, tmp$X)
        pressy[ioax, ioay] <- base::sum( base::sum((YScaleObjTest$matrix - 
                                                      modelPredy$Yhat)**2))
        pressyVars[ioax, ioay] <- base::sum((YScaleObjTest$matrix - 
                                               modelPredy$Yhat)**2)
        
        if (icv == 1) {
          pressyTot[ioax, ioay] <- pressy[ioax, ioay]
          pressyVarsTot[ioax, ioay] <- pressyVars[ioax, ioay]
        } else {
          pressyTot[ioax,ioay] <- pressyTot[ioax,ioay]+pressy[ioax,ioay]
          pressyVarsTot[ioax,ioay] <- pressyVarsTot[ioax,ioay]+pressyVars[ioax,ioay]
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
    }
    AllYhat <- rbind(AllYhat, AllYhatind)
  } # end icv
  
  if(verbose){
    base::cat("Cross-validation is complete.                              ")
    utils::flush.console()
  }
  
  KtrTr <- K
  modelMain <- list()
  modelMain$Model <- koplsModel(K = KtrTr, Y = Y, A = A, nox = oax, 
                                     preProcK = preProcK, preProcY = preProcY)
  modelMain$cv$Yhat <- base::matrix(data = base::unlist(Yhat), 
                                    ncol = length(Yhat[1,]))
  modelMain$cv$AllYhat <- AllYhat
  modelMain$cv$Tcv <- modelMain$cv$Yhat %*% modelMain$Model$Cp %*% modelMain$Model$Bt[[oax + 1]]
  modelMain$cv$Q2Yhat <- base::matrix(data = 0, nrow = oax+1, ncol = 1)
  modelMain$cv$Q2YhatVars <- base::matrix(data = 0, nrow = oax+1, ncol = 1)
  for(ioax in 1:(oax+1)){
    for(ioay in 1:1){
      modelMain$cv$Q2Yhat[ioax,ioay] <- 1 - pressyTot[ioax,ioay]/ssyTot
      modelMain$cv$Q2YhatVars[ioax,ioay] <- 1 - pressyVarsTot[ioax,ioay]/ssyVarsTot
    }
  }
  modelMain$cv$cvTestIndex <- cvTestIndex
  modelMain$cv$cvTrainingIndex <- cvTrainingIndex
  
  if (modelType == "da") {
    # Get sens/spec for each y-orth component
    for (i in 1:(oax + 1)) {
      if (drRule == "max") {
        predClass <- koplsMaxClassify(data = YhatDaSave[[i]])
      } else if (drRule == "fixed") {
        predClass <- koplsBasicClassify(data = YhatDaSave[[i]], 
                                        k = 1/nclasses)
      } else {
        warning(paste0('Decision rule given: ', drRule, 
                      ' is not valid/implemented.'))
      }
      
      # Calculate sensitivity and specificity
      daMetrics <- koplsSensSpec(trueClass = base::matrix(classVect[cvTestIndex]), 
                                 predClass = predClass)
      daMetrics$sens[i] <- daMetrics[["classResults"]][[i]][["sens"]]
      daMetrics$spec[i] <- daMetrics[["classResults"]][[i]][["spec"]]
      daMetrics$classvec[i] <- daMetrics[["classResults"]][[i]][["class"]]
      daMetrics$tot_sens[i] <- daMetrics[["totalResults"]][["sensTot"]]
      daMetrics$meanSens[i] <- daMetrics[["totalResults"]][["meanSens"]]
      daMetrics$meanSpec[i] <- daMetrics[["totalResults"]][["meanSpec"]]
    }
    
    # Calculate sensitivity and specificity
    daMetrics$confusMatrix <- koplsConfusionMatrix(true_class = classVect[cvTestIndex], 
                                                   pred = predClass)
    daMetrics$trueClass <- classVect[cvTestIndex]
    daMetrics$nclasses <- nclasses
    modelMain$da <- daMetrics
    modelMain$da$predClass <- predClass
    modelMain$da$decisionRule <- drRule
    
    # Change to original order if NFOLD CV
    if (cvType == "nfold") {
      cvOrder <- base::sort(x = cvTestIndex, decreasing = FALSE)
      modelMain$da$predClass <- modelMain$da$predClass[cvOrder]
      modelMain$da$trueClass <- modelMain$da$trueClass[cvOrder]
    }
  }
  
  # Change to original order if NFOLD CV
  if (cvType == "nfold") {
    cvOrder <- base::sort(x = cvTestIndex, decreasing = FALSE)
    modelMain$cv$Yhat <- modelMain$cv$Yhat[cvOrder, ]
    modelMain$cv$Tcv <- modelMain$cv$Tcv[cvOrder, ]
  }
  
  modelMain$class <- "koplscv"
  modelMain$da$args$oax <- oax
  modelMain$da$args$A <- A
  
  return(modelMain)
}


