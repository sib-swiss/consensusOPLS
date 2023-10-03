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
#' determine the optimal number of model components.
#' OUTPUT
% modelMain = Object with 'A' predictive components and 'oax'
%   Y-orthogonal components. Contains the following entries:
  %  cv = Cross-validation results:
    %  	 Q2Yhat = Total Q-square result for all Y-orthogonal components.
    %	 Q2YhatVars = Q-square result per Y-variable for all Y-orthogonal
    %       components.
    %	 Yhat = All predicted Y values as a concatenated matrix.
    %	 Tcv = Predictive score vector T for all cross-validation rounds.
    %	 cvTrainIndex = Indices for the training set observations during
    %       the cross-validation rounds.
    %	 cvTestIndex = Indices for the test set observations during the
    %       cross-validation rounds.
    %  da = Cross-validation results specifically for discriminant
    %       analysis (DA) cases:
      %    predClass = Predicted class list per class and Y-orthogonal
      %       components (integer values).
      %	 trueClass = Predicted class list per class and Y-orthogonal
      %       components (integer values).
      %	 sensSpec = Sensitivity and specificity values per class and
      %       Y-orthogonal components (integer values).
      %	 confusionMatrix = Confusion matrix during cross-validation
      %       rounds.
      %	 nclasses = Number of classes in model.
      %	 decisionRule = Decision rule used: 'max' or 'fixed'.
      %  args = Arguments to the function:
        %  	 A = Number of Y-predictive components.
        %	 oax = Number of Y-orthogonal components.
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
  if(!is.integer(A)){stop("A is not an integer.")}
  if(!is.integer(oax)){stop("oax is not an integer.")}
  if(!is.integer(nbrcv)){stop("nbrcv is not an integer.")}
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
  if(!is.integer(cvFrac)){stop("cvFrac is not an integer.")}
  if(!is.character(modelType){stop("modelType is not a character.")
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
    if(all(temp == c(0, 1))){
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
  # Yhat <- base::matrix(NA, nrow = nrow(Y), ncol = ncol(Y))
  # YhatDaSave=cell(1);
  # pressyVars=cell(1,1);
  # pressyVarsTot=cell(1,1);
  
  if(verbose){
    utils::txtProgressBar(min = 0, max = nbrcv, style = 3, 
                          width = 10, char = "=",
                          title = paste0('Please wait... cv round: 1 of ', nrcv))
  }
  
  AllYhat <- c()
  
  for(icv in 1:nbrcv){
    # Update progression bar
    if(verbose){
      print(paste0('Cross-validation round: ', icv,'...'))
      # Create a waitbar (you'll need the tcltk package)
      if (requireNamespace("tcltk", quietly = TRUE)) {
        h <- tcltk::tkProgressBar(title = 'Please wait...', 
                                  label = paste0('cv round: ', icv, 
                                                 ' of ', nrcv))
        tcltk::tkSetProgress(icv/nrcv*100, h)
      }
    }
    
    # Set up Cross-Validation
    cvSet <- koplsCrossValSet(K, Y, cvFrac, cvType, nbrcv, icv)
    cvTestIndex <- cvSet$testIndex
    cvTrainingIndex <- cvSet$trainingIndex
    
    # Get Kernel matrices 
    # % change so that this is done in the K matrix only once and 
    # selected by indices.
    KtrTr <- cvSet$KTrTr
    KteTe <- cvSet$KTeTe
    KteTr <- cvSet$KTeTr
    
    # Center Y and kernel matrices
    YScaleObj <- koplsScale(cvSet$yTraining,YcenterType,YscaleType)
    # I don't understand the useful of koplsScaleApply...
    YScaleObjTest <- koplsScale(cvSet$yTest,YScaleObj)
    
    if(preProcK == "mc"){
      KteTe <- koplsCenterKTeTe(KteTe,KteTr,KtrTr)
      KteTr <- koplsCenterKTeTr(KteTr,KtrTr)
      KtrTr <- koplsCenterKTrTr(KtrTr)
    }
    
    # Estimate K-OPLS model
    model <- koplsModel(KtrTr, YScaleObj$X, A, oax, 'no', 'no')
    
    # Set up model stats
    ssy <- sum( sum(YScaleObjTest$matrix)**2 )
    ssyVars <- sum(YScaleObjTest$matrix)**2
    ssx <- sum(diag(KteTe))
    
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
      for(ioay == 1){
        # Consensus OPLS predict Yhat
        modelPredy <- koplsPredict(KteTr, KteTe, KtrTr, model, ioax-1, 0)
        tmp <- koplsRescale(YScaleObj, modelPredy$Yhat)
        AllYhatind <- cbind(AllYhatind, tmp$X)
        pressy[ioax, ioay] <- sum(sum(YScaleObjTest$X - modelPredy$Yhat)**2)
        pressyVars[ioax, ioay] <- sum((YScaleObjTest$X - modelPredy$Yhat)**2)
        
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
            YhatDaSave[ioax, ioay] <- c()
          }
          
          # + mean on Yhat
          tmp <- koplsRescale(YScaleObj, modelPredy$Yhat)
          YhatDaSave[ioax,ioay] <- rbind(YhatDaSave[ioax, ioay], tmp$X)
        }
        
        # If highest number of oscs, save Yhat and Xhat
        if (ioax == oax+1) {  # && ioay == oay + 1) 
          if (icv == 1) {
            Yhat <- c()
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
      tcltk::tkSetProgress(icv, h)
    }
  }
  
  KtrTr <- K
  modelMain$koplsModel <- koplsModel(KtrTr, Y, A, oax, preProcK, preProcY)
  modelMain$cv$Yhat <- Yhat
  modelMain$cv$AllYhat <- AllYhat
  
  #is this correct?
  modelMain$cv$Tcv <- Yhat %*% modelMain$koplsModel$Cp %*% modelMain$koplsModel$Bt[oax + 1]
  
  modelMain.cv.Q2Yhat=[]; 
  modelMain.cv.Q2YhatVars=cell(1,1);
  
  
  for( ioax = 1:oax+1){
    for(ioay == 1){
      modelMain$cv$Q2Yhat(ioax,ioay) <- 1-pressyTot(ioax,ioay)/ssyTot
      modelMain$cv$Q2YhatVars(ioax,ioay) <- 1-pressyVarsTot(ioax,ioay)/ssyVarsTot
    }
  }
  
  
  # ---------------------------------------------------------
  # ----- MATLAB code

  #[scaleY]=koplsScale(Y,YcenterType,YscaleType);
  #KtrTr=koplsKernel(X,X,kernelType,kernelParam);
  #if (strcmp(preProcK,'mc'))
  #    KtrTr=koplsCenterKTrTr(KtrTr);
  #end
  #modelMain.koplsModel=koplsModel(KtrTr,scaleY.X,A,oax,preProcK,preProcY);


    modelMain.cv.cvTestIndex=cvTestIndex;
    modelMain.cv.cvTrainingIndex=cvTrainingIndex;

    if(strcmp(modelType,'da'))

        %get sens/spec for each y-orth component... eval of model
        for( i = 1:oax+1) %we would have no osc comps for dummy matrix...
                if(strcmp(drRule,'max'))
                    predClass=koplsMaxClassify(YhatDaSave{i,1});
                elseif(strcmp(drRule,'fixed'))
                    predClass=koplsBasicClassify(YhatDaSave{i,1},1/nclasses);
                else
                    warning(['Decision rule given: ',drRule,' is not valid/implemnted'])
                end
                %keyboard;
                [da.sensAllOsc{i}, da.specAllOsc{i}, da.classvecAllOsc{i}, da.tot_sensAllOsc{i},da.meanSensAllOsc{i},da.meanSpecAllOsc{i}]=koplsSensSpec(classVect(cvTestIndex), predClass);
        end
        

        
        % get sens/spec for max number of oscs.... (hmm redundant).
        
        if(strcmp(drRule,'max'))
            predClass=koplsMaxClassify(Yhat);
        elseif(strcmp(drRule,'fixed'))
            predClass=koplsBasicClassify(Yhat,1/nclasses);
        else
            warning(['Decision rule given: ',drRule,' is not valid/implemnted'])
        end
        

           [da.sens, da.spec, da.classvec, da.tot_sens,da.meanSens,da.meanSpec]=koplsSensSpec(classVect(cvTestIndex), predClass);
           [da.confusionMatrix]=koplsConfusionMatrix(classVect(cvTestIndex), predClass);
           da.trueClass=classVect(cvTestIndex);
            da.nclasses=nclasses;
        modelMain.da=da;
        modelMain.da.predClass=predClass;        
        modelMain.da.decisionRule=drRule;
                %CHANGE TO ORIGNAL ORDER IF NFOLD CV - for backward
                %compatibility and comparison w/ simca-p etc
        if(strcmp(cvType,'nfold'))
            [tmp,cvOrder]=sort(cvTestIndex);
            modelMain.da.predClass=modelMain.da.predClass(cvOrder);
            modelMain.da.trueClass=modelMain.da.trueClass(cvOrder);           
        end
        
    end


    %CHANGE TO ORIGNAL ORDER IF NFOLD CV - for backward
    %compatibility and comparison w/ simca-p etc
    if(strcmp(cvType,'nfold'))
        [tmp,cvOrder]=sort(cvTestIndex);

           modelMain.cv.Yhat=modelMain.cv.Yhat(cvOrder,:);
           modelMain.cv.Tcv=modelMain.cv.Tcv(cvOrder,:);

    end

if (verbose)
    close(h);
end

modelMain.release=release;
modelMain.args.oax=oax;
%modelMain.args.oay=oay;
modelMain.args.A=A;
modelMain.class='koplscv';

end






