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
#' cross-validation (integer). Only applicable for `mccv` or `mccvb` 
#' cross-validation (see `cvType`). 
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
#' `cv`: list. The cross-validation results.
#'      `Q2Yhat`: Total Q-square result for all Y-orthogonal components.
#'      `Q2YhatVars`: Q-square result per Y-variable for all Y-orthogonal
#'      components.
#'      `Yhat`: All predicted Y values as a concatenated matrix.
#'      `Tcv`: Predictive score vector T for all cross-validation rounds.
#'      `cvTrainIndex`: Indices for the training set observations during the 
#'      cross-validation rounds.
#'      `cvTestIndex`: Indices for the test set observations during the
#'      cross-validation rounds.
#' `da`: list. Cross-validation results specifically for discriminant analysis:
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
  
  # ----- Import local functions used in this code
  #1. koplsDummy
  #2. koplsReDummy
  #3. koplsCrossValSet
  #4. koplsScale
  #5. koplsCenterKTeTe
  #6. koplsCenterKTeTr
  #7. koplsCenterKTrTr
  #8. koplsModel
  #9. koplsPredict 
  #10. koplsRescale
  #11. koplsMaxClassify
  #12. koplsBasicClassify
  #13. koplsConfusionMatrix
  #14. koplsSensSpec
  
  # ----- Default values control
  if(cvType == "mccvb" & modelType != "da"){
    stop("Class balanced MC CV only applicable to `da` modelling.")
  } 
  if(!is.na(cvFrac) & !(cvType %in% c("mccv", "mccvb"))){
    stop("cvFrac only applicable to `mccv` or `mccvb` cvType.")
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
  YhatDaSave <- base::matrix(data = list(), nrow = oax+1, ncol = 1)

  if(verbose){
    utils::txtProgressBar(min = 0, max = nbrcv, style = 3, 
                          width = 10, char = "=",
                          title = paste0('Please wait... cv round: 1 of ', nrcv))
    # cette progress bar doit etre changée car ouvre une fenetre a part..
    # il faudrait que tout soit intégré à R
  }
  
  AllYhat <- c()
  
  for(icv in 1:nbrcv){
    # Update progression bar
    if(verbose){
      print(paste0('Cross-validation round: ', icv,' of ', nbrcv, ' ...'))
      # Create a waitbar (you'll need the tcltk package)
      if (requireNamespace("tcltk", quietly = TRUE)) {
        h <- tcltk::tkProgressBar(title = 'Please wait...', 
                                  label = paste0('cv round: ', icv, 
                                                 ' of ', nrcv))
        tcltk::setTkProgressBar(pb = h, value = icv/nrcv*100)
      }
    }
    
    # Set up Cross-Validation
    cvSet <- koplsCrossValSet(K = K, Y = Y, modelFrac = cvFrac, type = cvType, 
                              nfold = nbrcv, nfoldRound = icv)
    cvTestIndex <- cvSet$testIndex
    cvTrainingIndex <- cvSet$trainingIndex
    
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
            YhatDaSave[[ioax, ioay]] <- base::matrix(data = NA, 
                                                     nrow = nrow(modelPredy$Yhat),
                                                     ncol = ncol(modelPredy$Yhat))
          }
          
          # + mean on Yhat
          tmp <- koplsRescale(YScaleObj, modelPredy$Yhat)
          YhatDaSave[[ioax, ioay]] <- tmp$X
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
  
  if (verbose) {
    if (requireNamespace("tcltk", quietly = TRUE)) {
      h <- tcltk::tkProgressBar(title = 'finishing up...', 
                                label = '', max = nrcv, initial = icv)
      tcltk::setTkProgressBar(pb = h, value = icv)
    }
  }
  
  KtrTr <- K
  modelMain <- list()
  modelMain$koplsModel <- koplsModel(K = KtrTr, Y = Y, A = A, nox = oax, 
                                     preProcK = preProcK, preProcY = preProcY)
  modelMain$cv$Yhat <- Yhat
  modelMain$cv$AllYhat <- AllYhat
  modelMain$cv$Tcv <- Yhat %*% modelMain$koplsModel$Cp %*% modelMain$koplsModel$Bt[[oax + 1]]
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
        predClass <- koplsMaxClassify(data = YhatDaSave[[i,1]])
      } else if (drRule == "fixed") {
        predClass <- koplsBasicClassify(data = YhatDaSave[[i,1]], 
                                        k = 1/nclasses)
      } else {
        warning(paste0('Decision rule given: ', drRule, 
                      ' is not valid/implemented.'))
      }
      
      # Calculate sensitivity and specificity
      ########## BUG ICI ##########
      daMetrics <- koplsSensSpec(trueClass = classVect[cvTestIndex], 
                                 predClass = predClass)
      da$sensAllOsc[i] <- daMetrics$sens
      da$specAllOsc[i] <- daMetrics$spec
      da$classvecAllOsc[i] <- daMetrics$classvec
      da$tot_sensAllOsc[i] <- daMetrics$tot_sens
      da$meanSensAllOsc[i] <- daMetrics$meanSens
      da$meanSpecAllOsc[i] <- daMetrics$meanSpec
    }
    
    # Get sens/spec for max number of oscs
    if (drRule == "max") {
      predClass <- koplsMaxClassify(Yhat)
    } else if (drRule == "fixed") {
      predClass <- koplsBasicClassify(Yhat, 1/nclasses)
    } else {
      warning(paste0('Decision rule given: ', drRule, 
                     ' is not valid/implemented.'))
    }
    
    # Calculate sensitivity and specificity
    daMetrics <- koplsSensSpec(classVect[cvTestIndex], predClass)
    da$sens <- daMetrics$sens
    da$spec <- daMetrics$spec
    da$classvec <- daMetrics$classvec
    da$tot_sens <- daMetrics$tot_sens
    da$meanSens <- daMetrics$meanSens
    da$meanSpec <- daMetrics$meanSpec
    da$confusionMatrix <- koplsConfusionMatrix(classVect[cvTestIndex], predClass)
    da$trueClass <- classVect[cvTestIndex]
    da$nclasses <- nclasses
    modelMain$da <- da
    modelMain$da$predClass <- predClass
    
    # Change to original order if NFOLD CV
    if (cvType == "nfold") {
      cvOrder <- order(cvTestIndex)
      modelMain$da$predClass <- modelMain$da$predClass[cvOrder]
      modelMain$da$trueClass <- modelMain$da$trueClass[cvOrder]
    }
  }
  
  # Change to original order if NFOLD CV
  if (cvType == "nfold") {
    cvOrder <- order(cvTestIndex)
    modelMain_cv_Yhat <- modelMain$cv$Yhat[cvOrder, ]
    modelMain_cv_Tcv <- modelMain_cv_Tcv[cvOrder, ]
  }
  
  if (verbose) {
    close(h)
  }
  
  return(modelMain = list("release" = release,
                          "cv" = list("Q2Yhat" = xx,
                                      "Q2YhatVars" = xx,
                                      "Yhat" = modelMain_cv_Yhat,
                                      "Tcv" = modelMain_cv_Tcv,
                                      "cvTrainIndex" = xx,
                                      "cvTestIndex" = xx),
                          "da" = list("predClass" = xx,
                                      "trueClass" = xx,
                                      "sensSpec" = xx,
                                      "confusionMatrix" = xx,
                                      "nclasses" = xx,
                                      "decisionRule" = xx,
                                      "arg" = list("oax" = oax,
                                                   "A" = A)),
                          "class" = "koplscv"))
}


