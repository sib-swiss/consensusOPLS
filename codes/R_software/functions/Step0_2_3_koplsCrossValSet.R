#' koplsCrossValSet
#' Generates set of training/test observations for the Cross-Validation
#'
#' @param K: Kernel matrix. 
#' @param Y: Response matrix.
#' @param modelFrac
#' @param type: Type of cross-validation. It can be `nfold` for n-fold, 
#' `mccv` for Monte Carlo CV, `mccvb` for Monte Carlo class-balanced CV.
#' Default is `mccv`.
#' @param nfold 
#' @param nfoldRound 
#'
#' @return
#' Object `cvSet` with the following entries:
#' `type` = Cross-validation type.
#' `nfold` = Number of nfold rounds.
#' `nfoldRound` = The current nfold round.
#' `KTrTr` = Kernel training matrix; KTrTr = <phi(Xtr),phi(Xtr)>.
#' `KTeTr` = Kernel test/training matrix; KTeTr = <phi(Xte),phi(Xtr)>.}
#' `KTeTe` = Kernel test matrix; KTeTe = <phi(Xte),phi(Xte)>.
#' `yTraining` = Y training set.
#' `yTest` = Y test set.
#' `trainingIndex` = Indices of training set observations.
#' `testIndex` = Indices of test set observations.
#'
#' @examples
#' TO DO

koplsCrossValSet <- function(K, Y, modelFrac, type, 
                             nfold = NA, nfoldRound = NA){
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
      warning("type is not nfold, nfold is not defined (missing argument).")
    }
  }
  if(!is.na(nfoldRound)){
    if(!is.numeric(nfoldRound)){stop("nfoldRound is not a number.")
    }
    if(type != "nfold"){
      nfoldRound <- NA
      warning("type is not nfold, nfoldRound is not defined (missing argument).")
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
    ##### ----- code that I do
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
}