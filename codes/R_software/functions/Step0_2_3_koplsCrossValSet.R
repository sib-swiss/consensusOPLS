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
#' Y <- matrix(1:100, nrow = 10)
#' K <- matrix(1:100, nrow = 10)  
#' type <- "nfold"  
#' modelFrac <- 0.7  
#' nfold <- 5  
#' nfoldRound <- 1 
#' test <- koplsCrossValSet(K = K, Y = Y, modelFrac = modelFrac,
#'                          type = type, nfold = nfold, nfoldRound = nfoldRound) 

koplsCrossValSet <- function(K, Y, modelFrac, type, 
                             nfold = NA, nfoldRound = NA){
  # ----- Variable format control
  if(!is.matrix(K)){stop("K is not a matrix.")}
  if(!is.matrix(Y)){stop("Y is not a matrix.")}
  if(!is.numeric(modelFrac)){stop("modelFrac is not numeric.")}
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
    if(all(temp == c(0, 1))){
      classVect <- koplsReDummy(Y)
    } else{
      classVect <- Y
    }
    
    # Find all class labels
    minset <- unique(classVect)
    
    # For each class
    for(i in 1:length(minset)){
      # Find samples of current class
      ind <- which(classVect == minset[i])
      
      # Randomize
      ran <- runif(length(ind))
      rand_ind <- order(ran)
      ind <- ind[rand_ind]
      
      # Adjusted model with randomized vector
      modelLim <- ceiling(length(ind) * modelFrac)
      modInd <- c(modInd, ind[1:modelLim])
      predInd <- c(predInd, ind[(modelLim + 1):length(ind)])
    }
  }
  
  # Define Monte-Carlos Cross Validation
  if(type == "mccv"){
    # Randomize
    ran <- runif(length(K[, 1]))
    rand_ind <- order(ran)
    ind <- 1:length(ran)
    ind <- ind[rand_ind]
    
    # Adjusted model with randomized vector
    modelLim <- ceiling(length(ind) * modelFrac)
    modInd <- ind[1:modelLim]
    predInd <- ind[(modelLim + 1):length(ind)]
  }
  
  # Define N-fold Cross Validation
  if(type == "nfold"){
    predInd <- base::seq(from = nfoldRound, to = nfold, by = nfold)
    modInd <- setdiff(1:length(Y[, 1]), predInd)
  }
  
  # Apply cross validation
  if(nrow(Y) == ncol(Y)){
    KTrTr <- K[modInd, modInd]
    KTeTr <- K[predInd, modInd]
    KTeTe <- K[predInd, predInd]
  } else{
    KTrTr <- NA
    KTeTr <- NA
    KTeTe <- NA
  }
  
  # Return the final CV Set
  return(cvSet = list("type" = type,
                      "nfold" = nfold,
                      "nfoldRound" = nfoldRound,
                      "KTrTr" = KTrTr,
                      "KTeTr" = KTeTr,
                      "KTeTe" = KTeTe,
                      "yTraining" = Y[modInd, ],
                      "yTest" = Y[predInd, ],
                      "trainingIndex" = modInd,
                      "testIndex" = predInd))
}
