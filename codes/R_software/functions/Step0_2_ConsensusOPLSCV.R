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
#' @export
#'
#' @examples
#' 

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
  
  # ----- Define local functions used in this code
  koplsDummy <- function(class, numClasses = NA){
    # class: integer vector with classes to define.
    # numClasses: pre-defined the number of classes in the output. By default
    # is equal to the length of `class` vector.
    
    # Define parameters
    labels <- base::unique(class)
    labels_sorted <- base::sort(labels, decreasing = FALSE)
    samples_number <- base::length(class)
    
    # Variable format control
    if(is.na(class)){stop("class must be contain an integer vector.")}
    if(is.na(numClasses)){
      sample_labels <- base::length(labels)
      ncol <- sample_labels
      default <- TRUE
    } else{
      sample_labels <- base::cut(x = class, breaks = numClasses,
                                 ordered_result = TRUE)
      ncol <- base::length(levels(sample_labels))
      default <- FALSE
    }
    
    # Matrix initialization
    dummy <- base::matrix(data = 0, nrow = samples_number, ncol = ncol)
    
    # Search for class membership
    for(i in 1:ncol){
      if(default){
        dummy[class %in% labels_sorted[i], i] <- 1
      } else{
        labels_cut <- cbind(lower <- as.numeric(base::sub("\\((.+),.*", "\\1", 
                                                          levels(sample_labels)[i],
                                                          fixed = FALSE)),
                            upper <- as.numeric(base::sub("[^,]*,([^]]*)\\]", "\\1", 
                                                          levels(sample_labels)[i],
                                                          fixed = FALSE)))
        dummy[(class > labels_cut[1, "lower"] & 
               class <= labels_cut[1, "upper"]), i] <- 1
      }
    }
    # Return a list with 2 elements
    # matrix: A matrix with rows corresponding to observations and columns 
    # to classes. Each element in matrix is either one (observation belongs to 
    # class) or zero (observation does not belong to class).
    # labels_sorted: The class labels that are found in class in sorted order.
    return(list("matrix" = dummy,
                "labels_sorted" = labels_sorted))
    # example:
    # class <- c(5, 1, 2, 3, 4, 3, 2, 4, 3, 1, 3)
    # numClasses <- 2
  }
  koplsReDummy <- function(Y){
    # Y: Dummy matrix.
    
    # Extract size of matrix
    n <- nrow(Y)
    m <- ncol(Y)
    
    # Create an empty vector
    classVect <- base::rep(x = NA, times = n)
    
    # Rebuild the vector
    for(i in 1:m){
      classVect[Y[, i] == 1] <- i
    }
    
    # Return a (integer) class vector from a binary matrix
    return(classVect)
    # example:
    # Y <- koplsDummy(class = c(5, 1, 2, 3, 4, 3, 2, 4, 3, 1, 3),
    #                 numClasses = NA)$matrix
  }
  koplsCrossValSet <- function(K, Y, modelFrac, type, 
                               nfold = NA, nfoldRound = NA){
    # Generates set of training/test observations for the CV
    # How the sets are generated is determined by the 'type' parameter, 
    # which can be either 'nfold' for n-fold cross-validation, 
    # 'mccv' for Monte Carlo CV, 'mccvb' for Monte Carlo class-balanced CV. 
    
    # ----- Variable format control
    if(!is.matrix(K)){stop("K is not a matrix.")}
    if(!is.matrix(Y)){stop("Y is not a matrix.")}
    #if(modelFrac){}
    if(!is.character(type)){stop("type is not a character.")
    } else{
      if(!(type %in% c("nfold", "mccv", "mccvb"))){
        stop("type must be `nfold`, `mccv` or `mccvb`.")
      }
    }
    if(!is.na(nfold)){
      if(!is.numeric(nfold)){stop("nfold is not a number.")
      }
      if(type != "nfold"){
        nfold <- NA
        warning("type is not nfold, nfold is defined as missing.")
      }
    }
    if(!is.na(nfoldRound)){
      if(!is.numeric(nfoldRound)){stop("nfoldRound is not a number.")
      }
      if(type != "nfold"){
        nfoldRound <- NA
        warning("type is not nfold, nfoldRound is defined as missing.")
      }
    }
    
    # Define Monte-Carlos Cross Validation - class Balanced
    if(type == "mccvb"){
      # check if Y is dummy or labels
      temp <- unique(Y)
      if(temp %in% c(0, 1)){
        classVect <- koplsReDummy(Y)
      } else{
        classVect <- Y
      }
      
      # Stock class labels
      minset <- unique(ClassVect)
      for(i in 1:length(minset)){
        # Find samples of current class
        ind <- (classVect == minset[i])
        
        # Randomization
        ##### ----- code que je ne sais pas traduire
        # %randomize
        # ran=rand(length(ind),1); %randomize
        # [tmp,rand_ind]=sort(ran); %sort randomize number to get randomized index string
        # ind=ind(rand_ind); %apply randomization on the real index vector
        # %-end randomize
        ##### ----- fin
        rand <- sample(ind)
        
        # Adjusted model with randomized vector
        modelLim <- round(length(ind)*modelFrac)
        modInd <- ind[1:modelLim]
        predInd <- ind[(modelLim+1):length(classVect)]
      }
    }
    
    # Define Monte-Carlos Cross Validation
    if(type == "mccv"){
      # Randomization
      ##### ----- code que je ne sais pas traduire
      # %randomize
      # ran=rand(length(K(:,1)),1); %randomize
      # [tmp,rand_ind]=sort(ran); %sort randomize number to get randomized index string
      # ind=[1:length(ran)]';
      # ind=ind(rand_ind); %apply randomization on the real index vector
      # %-end randomize
      ##### ----- fin
      rand <- sample(1:nrow(K))
      
      # Adjusted model with randomized vector
      modelLim <- round(length(ind)*modelFrac)
      modInd <- ind[1:modelLim]
      predInd <- ind[(modelLim+1):nrow(K)]
    }
    
    # Define N-fold Cross Validation
    if(type == "nfold"){
      predInd <- nfoldRound:nfold:nrow(Y)
      modInd <- setdiff(1:nrow(Y), predInd)
    }
    
    # Apply cross validation
    if(nrow(Y) == ncol(Y)){
      cvSet.KTrTr <- K(modInd,modInd)
      cvSet.KTeTr <- K(predInd,modInd)
      cvSet.KTeTe <- K(predInd,predInd)
    } else{
      cvSet.KTrTr <- NA
      cvSet.KTeTr <- NA
      cvSet.KTeTe <- NA
    }
    
    # Return the final CV Set
    return(cvSet = list("type" = type,
                        "nfold" = nfold,
                        "nfoldRound" = nfoldRound,
                        "yTraining" = Y[modInd, ],
                        "yTest" = Y[predInd, ],
                        "trainingIndex" = modInd,
                        "testIndex" = predInd))
    # example:
    # TO DO
  }
  koplsScale <- function(X, centerType, scaleType){
    # Mean-centering and scaling matrix
    # X the matrix to be centered and scaled
    # centerType for mean (or no) centering
    # scaleType for unit variance scaling, Pareto scaling of no scaling
    
    # Variable format control
    if(!is.matrix(X)){stop("X is not a matrix.")}
    if(!is.character(centerType)){stop("centerType is not a character.")
    } else{ 
      if(!(centerType %in% c("mc", "no"))){
        stop("centerType must be `mc` or `no`.")}
    }
    if(!is.character(scaleType)){stop("scaleType is not a character.")
    } else{ 
      if(!(scaleType %in% c("uv", "pa", "no"))){
        stop("scaleType must be `uv`, `pa` or `no`.")}
    }
    
    # Define variables
    m <- nrow(X)
    
    # Center the matrix
    if(centerType == "mc"){
      X <- apply(X, 2, FUN = function(X){X - mean(X)})
    }
    
    # Scale the matrix
    if(scaleType == "uv"){
      X <- apply(X, 2, FUN = function(X){X/sd(X)})
    }
    if(scaleType == "pa"){
      X <- apply(X, 2, FUN = function(X){X/sqrt(sd(X))})
    }
    
    # Return 
    return(list("centerType" = centerType,
                "scaleType" = scaleType,
                "meanV" = mean(X),
                "stdV" = sd(X),
                "matrix" = X))
    # example:
    # X <- matrix(c(1,4,7, 8,4,0, 3,6,9), nrow=3)
  }
  koplsCenterKTeTe <- function(KteTe, KteTr, KtrTr){
    # Centering function for test Kernel
    # KteTe: Test kernel matrix, KteTe = <phi(Xte), phi(Xte)>
    # KteTr: Test/training kernel matrix, KteTr = <phi(Xte), phi(Xtr)>
    # KtrTr: raining kernel matrix, KtrTr = <phi(Xtr), phi(Xtr)>
    
    # Variable format control
    if(!is.matrix(KteTe)){stop("KteTe is not a matrix.")}
    if(!is.matrix(KteTr)){stop("KteTr is not a matrix.")}
    if(!is.matrix(KtrTr)){stop("KtrTr is not a matrix.")}
    
    # Define parameters
    Itrain <- base::diag(ncol(KteTr));
    I_nTrain <- matrix(1, nrow = ncol(KteTr), ncol = 1)
    nTrain <- ncol(KteTr);
    
    I <- base::diag(nrow(KteTr));
    I_n <- matrix(1, nrow = nrow(KteTr), ncol = 1)
    n <- nrow(KteTr);
    
    # Center the kernel
    D_te = (1/nTrain)%*%I_n%*%I_nTrain;
    KteTe = KteTe - D_te%*%KteTr - KteTr%*%D_te + D_te%*%KtrTr%*%D_te;
    
    # Return the centered test kernel matrix
    return(KteTe)
    # example:
    # 
  }
  koplsCenterKTeTr <- function(KteTr, KteTr){
    # Centering function for the hybrid test/training kernel, which
    # is constructed from the test matrix Xte and the training matrix
    # Xtr as KteTr = <phi(Xte), phi(Xtr)>.

    # Variable format control
    if(!is.matrix(KteTr)){stop("KteTr is not a matrix.")}
    if(!is.matrix(KtrTr)){stop("KtrTr is not a matrix.")}
    
    # Define parameters
    Itrain <-  base::diag(nrow(KtrTr))
    I_nTrain <- matrix(1, nrow = nrow(KtrTr), ncol = 1)
    nTrain <- nrow(KtrTr)
    
    I <- base::diag(nrow(KteTr))
    I_n <- matrix(1, nrow = nrow(KteTr), ncol = 1)
    n <- nrow(KteTr)
    
    KteTr = (KteTr-(1/nTrain)%*%I_n%*%I_nTrain %*% KtrTr) %*% 
      (Itrain-(1/nTrain)%*%I_nTrain%*%I_nTrain)
    
    # Return
    return(KteTr)
    # example:
    #
  }
  koplsCenterKTrTr <- function(K){
    # Centering function for the training kernel, which is constructed
    # from the training matrix Xtr as K = <phi(Xtr), phi(Xtr)>
    # K: training kernel matrix; K = <phi(Xtr), phi(Xtr)>.
    I <- base::diag(nrow(K))
    I_n <- matrix(1, nrow = nrow(K), ncol = 1)
    K = (I- (1/nrow(K))%*% I_n%*%I_n) %*% K%*%(I-(1/n)%*%I_n%*%I_n);
    
    # Return the centered kernel matrix
    return(K)
    # example:
    #
  }
  koplsModel <- function(K, Y, A, nox, preProcK, preProcY){
    # Function for training a K-OPLS model. The function constructs a
    # predictive regression model for predicting the values of 'Y by
    # using the information in 'K'. The explained variation is separated
    # into predictive components (dimensionality is determined by the
    # parameter 'A') and 'Y'-orthogonal components (dimensionality
    # determined by the parameter 'nox').
    
    # K: Kernel matrix (un-centered); K = <phi(Xtr),phi(Xtr)>.
    # Y: Response matrix (un-centered/scaled).
    # A: Number of predictive components.
    # nox: Number of Y-orthogonal components.
    # preProcK: Pre-processing parameters for the 'K' matrix:
    #     'mc' for mean-centering, 'no' for no centering.
    # preProcY = Pre-processing parameters for the 'Y' matrix:
    #     'mc' for mean-centering, 'uv' for mc + scaling to unit variance,
    #     'pa' for mc + Pareto, 'no' for no scaling.

    # Variable format control
    if(!is.matrix(K)){stop("K is not a matrix.")}
    if(!is.matrix(Y)){stop("Y is not a matrix.")}
    if(!is.numeric(A)){stop("A is not a numeric.")}
    if(!is.numeric(nox)){stop("nox is not a numeric.")}
    if(!is.character(preProcK)){stop("preProcK is not a character.")
    } else{
      if(!(preProcK %in% c("mc", "no"))){stop("preProcK must be `mc` or `no`.")}
    }
    if(!is.character(preProcY)){stop("preProcY is not a character.")
    } else{ 
      if(!(preProcY %in% c("mc", "uv", "pareto", "no"))){
        stop("preProcY must be `mc`, `uv`, `pareto` or `no`.")}
    }
    
    # Initialise parameters
    I <- diag(nrow = norw(K))
    Kmc <- K
    
    # Check kernel centering
    if(preProcK == "mc"){
      Kmc <- koplsCenterKTrTr(K)
    }
    K <- matrix(list(), nrow = nox+1, ncol = nox+1)
    K[1,1] <- Kmc
    
    # Check response matrix
    Y_old <- Y
    if(preProcY != "no"){
      if(preProcY == "mc"){
        scale <- "no"
      } else{
        scale <- preProcY
      }
      scaleParams <- koplsScale(Y, "mc", scale)
      Y <- scaleParams[[X]]
    }
    
    # KOPLS model estmation
    ## step 1
    svd_obj <- base::svd(t(Y)%*%K{1,1}%*%Y)
    Cp <- svd_obj$d[, 1:A]
    Sp <- svd_obj$u[1:A, 1:A]
    
    ## step 2
    Up <- Y%*%Cp
    
    ## step3
    Tp <- c() ; Bt <- c() ; co <- c()
    to <- c() ; toNorm <- c()
    for(i in 1:nox){
      ## step 4
      Tp <- base::append(t(K[1, i])%*%Up%*%Sp**{-1/2})
      Bt <- append( solve( t(Tp[i])%*%Tp[i]) %*%t(Tp[i])%*%Up)
      
      ## step 5
      temp <- base::svd( t(Tp[i])%*% (K[i,i]-Tp[i]%*%t(Tp[i])) %*% Tp[i]) 
      co <- append(temp$d[, 1])
      so <- append(temp$u[1,1])
      
      ## step 6
      to <- append( (K[i,i]-Tp[i]%*%t(Tp[i])) %*%Tp[i]%*%co[i]%*%(so[i]**{-1/2}) )
      
      ## step 7
      toNorm <- append( sqrt(t(to[i])%*%to[i]) )
      
      ## step 8
      to[i] <- to[i]/toNorm[i]
      
      ## step 9
      K[1, i+1] <- K[1,i]%*%(I-to[i]%*%t(to[i]))
      
      ## step 10
      K[i+1, i+1] <- (I-to[i]%*%t(to[i]))%*%K[i, i]%*%(I-to[i]%*%t(to[i]))
    } ## step 11
    
    ## step 12
    Tp[nox+1] <- t(K[1, nox+1])%*%Up%*%(Sp**{-1/2})
    
    ## step 13
    Bt[i+1] <- solve( t(Tp[nox+1])%*%Tp[nox+1] )%*%t(Tp[nox+1])%*%Up
    
    # ---------- extra stuff -----------------
    # should work but not fully tested (MB 2007-02-19)
    sstot_Y = sum(sum(Y.*Y));
    F=Y-Up*t(Cp)
    R2Y = 1 - sum(sum( F.*F))/sstot_Y
    # --------- #
    
    EEprime <- K[nox+1, nox+1]-Tp[nox+1]%*%t(Tp[nox+1])
    sstot_K <- sum(diag(K[1,1]))
    
    R2X <- c() ; R2X0 <- c() ; R2XC <- c() ; R2Yhat <- c()
    for(i in 1:(nox+1)){
      rss <- sum( diag(K[i,i]-Tp[i]%*%t(Tp[i])) )
      R2X <- append(1- rss/sstot_K)
      
      rssc <- sum( diag(K[1,1]-Tp[i]%*%t(Tp[i])) )
      R2XC <- append(1- rssc/sstot_K)
      
      # R2Yhat 22 Jan 2010 / MR - not fully tested
      Yhat <- Tp[i] %*% Bt[i] %*% t(Cp)
      R2Yhat <- append(1 - sum(sum((Yhat-Y)**2))/sstot_Y] )
    } # fin K-OPLS model
    
    # Convert to matrix structure
    To <- matrix(0, nrow = nrow(Tp), ncol = nox)
    for(i in 1:length(to)){
      To[, i]<- to[i]
    }
    
    return(model = list("Cp" = Cp, "Sp" = Sp, "Up" = Up,
                        "Tp" = Tp, "T" = Tp[nox+1], "co" = co,
                        "so" = so, "to" = to, "toNorm" = toNorm,
                        "Bt" = Bt, "A" = A, "nox" = nox, "K" = K,
                        "EEprime" =EEprime, "sstot_K" = sstot_K,
                        "R2X" = R2X, "R2XO" = R2XO, "R2XC" = R2XC,
                        "sstot_Y" = sstot_Y, "R2Y" = R2Y,
                        "R2Yhat" = R2Yhat, # R2Yhat 22 Jan 2010 / MR
                        "preProcK" = preProcK, "preProcY" = preProcY,
                        "preProcParamsY" = scaleParams))

  }
  
  
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
    if(temp %in% c(0, 1)){
      if(ncol(Y) == 1){Y <- koplsDummy(Y)}
      classVect <- koplsReDummy(Y)
    } elseif((temp %% 1) == 0) & ncol(Y) == 1){
      classVect <- Y
      Y <- koplsDummy(Y+1)
    } else{
      stop("modelType is `da`, but Y appears to be neither dummy matrix nor 
      a vector of (integer) class labels")
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
      utils::txtProgressBar(min = (icv-1)/nbrcv, max = nbrcv, style = 3, 
                            width = 10, char = "=",
                            title = paste0('Please wait... cv round: ',
                                           icv, ' of ', nrcv))
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
    model <- koplsModel(KtrTr,YScaleObj.X,A,oax,'no','no');
    
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
        
      }
    }
  }
  
  
  # ---------------------------------------------------------
  # ----- MATLAB code

%%

for( icv = 1:nbrcv)
	for( ioax = 1:oax+1)
        for( ioay = 1:1)%oay+1)        
            % ConsensuOPLS predict yYhat                                   
            [modelPredy]=koplsPredict(KteTr,KteTe,KtrTr, model,ioax-1,0);
            tmp=koplsRescale(YScaleObj,modelPredy.Yhat); 
            AllYhatind=[AllYhatind tmp.X];
            pressy(ioax,ioay)=sum(sum((YScaleObjTest.X-modelPredy.Yhat).^2));        
            pressyVars{ioax,ioay}=(sum((YScaleObjTest.X-modelPredy.Yhat).^2));	            
            
            if((icv==1))
                pressyTot(ioax,ioay)=pressy(ioax,ioay);
                pressyVarsTot{ioax,ioay}=pressyVars{ioax,ioay};
            
            else
                pressyTot(ioax,ioay)=pressyTot(ioax,ioay)+pressy(ioax,ioay);
                pressyVarsTot{ioax,ioay}=pressyVarsTot{ioax,ioay}+pressyVars{ioax,ioay};
             
            end
            
            
            %if 'da' save Yhat for all rounds
            if(strcmp(modelType,'da'))
                    if(icv==1)
                        YhatDaSave{ioax,ioay}=[];                        
                    end                    
                 
                    
                    %+mean on Yhat
                    tmp=koplsRescale(YScaleObj,modelPredy.Yhat);
                    YhatDaSave{ioax,ioay}=[YhatDaSave{ioax,ioay};tmp.X];                               
            end
            
            %if highest number of oscs - save Yhat and Xhat
            if(ioax==oax+1)% && ioay==oay+1)                
                    if(icv==1)
                        Yhat=[];
                      
                    end
					              
                       tmp=koplsRescale(YScaleObj,modelPredy.Yhat);
                       Yhat=[Yhat;tmp.X];
				
            end
            
        end
	
    end
    AllYhat=[AllYhat; AllYhatind];
	
end %end icv

if (verbose)
    waitbar(icv/(nbrcv),h,'finishing up...');
end

%[scaleY]=koplsScale(Y,YcenterType,YscaleType);
%KtrTr=koplsKernel(X,X,kernelType,kernelParam);
%if (strcmp(preProcK,'mc'))
%    KtrTr=koplsCenterKTrTr(KtrTr);
%end
%modelMain.koplsModel=koplsModel(KtrTr,scaleY.X,A,oax,preProcK,preProcY);

KtrTr=K;
modelMain.koplsModel=koplsModel(KtrTr,Y,A,oax,preProcK,preProcY);

modelMain.cv.Yhat=Yhat;
modelMain.cv.AllYhat=AllYhat;

%is this correct?
modelMain.cv.Tcv=Yhat*modelMain.koplsModel.Cp*modelMain.koplsModel.Bt{oax+1};

    modelMain.cv.Q2Yhat=[]; 
    modelMain.cv.Q2YhatVars=cell(1,1);


    for( ioax = 1:oax+1)
        for( ioay = 1:1)%oay+1)   
            modelMain.cv.Q2Yhat(ioax,ioay)=1-pressyTot(ioax,ioay)./ssyTot;
            modelMain.cv.Q2YhatVars{ioax,ioay}=1-pressyVarsTot{ioax,ioay}./ssyVarsTot;
        end
    end

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






