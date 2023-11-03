#' RVConsensusOPLSPerm
#' Y random permutations for model validity estimation, with leave-on-out CV.
#' 
#' @param data: list. The collection list containing each block of data.
#' @param Y: matrix. The response matrix.
#' @param nbruns: numeric. Number of random permutations. 
#' @param PredLVs: numeric. Number of predictive component(s).
#' @param maxOrtholvs: numeric. Maximum number of orthogonal LVs to compute.
#'
#' @return 
#' A list with the two results plots:
#' `plot_R2val`: plot. Plot RV vs R2val.
#' `plot_Q2val`: plot. Plot RV vs Q2val.
#'
#' @examples
#' TO DO

RVConsensusOPLSPerm <- function(data, Y, nbruns, PredLVs, 
                                maxOrtholvs){
  # Variable format control
  if(!is.list(data)){stop("data is not a list.")}
  if(!is.matrix(Y)){stop("Y is not a matrix.")}
  if(!is.numeric(nbruns)){stop("A is not numeric.")}
  if(!is.numeric(PredLVs)){stop("PredLVs is not numeric.")}
  if(!is.numeric(maxOrtholvs)){stop("maxOrtholvs is not numeric.")}

  # Function loading control
  if (!exists("RVConsensusOPLS", mode = "function")) {
    warning("Remember to load the source code for the `RVConsensusOPLS` function.")
  }
  if (!exists("RVmodified", mode = "function")) {
    warning("Remember to load the source code for the `RVmodified` function.")
  }
  
  # Init parameters
  NbObs <- nrow(data[[1]])
  PermRes <- list()
  
  # Extract parameters RVConsensusOPLS to define permutation parameters
  modelCV <- RVConsensusOPLS(data = collection, Y = Y, A = PredLVs, 
                             maxOrtholvs = maxOrtholvs, nrcv = NbObs,
                             cvType = "nfold", modelType = "da", 
                             verbose = FALSE)
  PermRes$lvnum[1] <- modelCV$model$cv$OrthoLVsOptimalNum + PredLVs
  PermRes$R2val[1] <- modelCV$model$Model$R2Yhat[base::length(modelCV$model$Model$R2Yhat)]
  PermRes$DQ2val[1] <- modelCV$model$cv$DQ2Yhat[PermRes$lvnum]
  PermRes$Q2val[1] <- modelCV$model$cv$Q2Yhat[PermRes$lvnum]
  PermRes$PredAc[1] <- modelCV$model$da$tot_sens[2]
  PermRes$Y[[1]] <- Y
  PermRes$RV[1] <- RVmodified(X = Y, Y = Y)
  
  # Progression bar
  base::cat("Please wait... The random permutation process begins.")
  utils::flush.console()
  
  # Permutations
  for(i in 1:nbruns){
    #Fix and reproduce the random
    set.seed(i)
    
    # Progress bar update
    cat("\r", "                                                         ", "\r")
    progress <- round(i * 100 / nbruns, 0)
    base::cat(base::sprintf("Progression : %.2f%% \r", progress))
    utils::flush.console()
    
    # Random permutation of Y rows
    Ys <- Y[base::sample(x = 1:nrow(Y), size = nrow(Y), replace = FALSE, prob = NULL), ]
    
    # Redo the Consensus OPLS-DA with RV coefficients weighting
    modelCV <- RVConsensusOPLS(data = collection, Y = Ys, A = PredLVs, 
                               maxOrtholvs = maxOrtholvs, nrcv = NbObs,
                               cvType = "nfold", modelType = "da", 
                               verbose = FALSE)
    PermRes$lvnum[i+1] <- modelCV$model$cv$OrthoLVsOptimalNum + PredLVs
    PermRes$R2val[i+1] <- modelCV$model$Model$R2Yhat[length(modelCV$model$Model$R2Yhat)]
    PermRes$DQ2val[i+1] <- modelCV$model$cv$DQ2Yhat[PermRes$lvnum[i+1]]
    PermRes$Q2val[i+1] <- modelCV$model$cv$Q2Yhat[PermRes$lvnum[i+1]]
    PermRes$PredAc[i+1] <- modelCV$model$da$tot_sens[2]
    PermRes$Y[[i+1]] <- Ys
    PermRes$RV[i+1] <- RVmodified(X = Y, Y = Ys)
  }
  
  # Progress bar update
  base::cat("Random permutation is complete.                                  ")
  utils::flush.console()
  
  # Plot the results
  PermRes_df <- data.frame(
    "RV" = PermRes$RV,
    "Values" = c(PermRes$R2val, PermRes$Q2val),
    "Type" = c(rep("RV", times = length(PermRes$R2val)),
               rep("Q2", times = length(PermRes$Q2val)))
  )
  
  # Scatterplots
  p1 <- ggplot2::ggplot(data = PermRes_df[which(PermRes_df$Type == "RV"), ],
                        aes(x = RV, y = Values)) +
    ggplot2::geom_point(size = 2.5)+
    ggplot2::geom_smooth(method = 'lm')+
    ggplot2::xlab("RV") +
    ggplot2::ylab("R2 values") +
    theme_graphs
  
  p2 <- ggplot2::ggplot(data = PermRes_df[which(PermRes_df$Type == "Q2"), ],
                        aes(x = RV, y = Values)) +
    ggplot2::geom_point(size = 2.5)+
    ggplot2::geom_smooth(method = 'lm')+
    ggplot2::xlab("RV") +
    ggplot2::ylab("Q2 values") +
    theme_graphs
  
  # Scatterplots all in one
  p3 <- ggplot2::ggplot(data = PermRes_df,
                        aes(x = RV, y = Values, col = Type)) +
    ggplot2::geom_point(size = 2.5)+
    ggplot2::xlab("RV") +
    ggplot2::ylab("R2 and Q2 values") +
    theme_graphs
  
  # Histograms
  p4 <- ggplot2::ggplot(data = PermRes_df[which(PermRes_df$Type == "RV"), ],
                        aes(x = Values)) +
    ggplot2::geom_histogram(color="black", fill="white")+
    ggplot2::geom_density(alpha=.2)+
    ggplot2::xlab("Frequency") +
    ggplot2::ylab("R2 values") +
    theme_graphs
  
  p5 <- ggplot2::ggplot(data = PermRes_df[which(PermRes_df$Type == "Q2"), ],
                        aes(x = Values)) +
    ggplot2::geom_histogram(color="black", fill="white")+
    ggplot2::geom_density(alpha=.2)+
    ggplot2::xlab("Frequency") +
    ggplot2::ylab("Q2 values") +
    theme_graphs
  
  # Histograms all in one
  p6 <- ggplot2::ggplot(data = PermRes_df,
                        aes(x = Values, col = Type)) +
    ggplot2::geom_histogram(size = 2.5)+
    ggplot2::geom_density(alpha=.2)+
    ggplot2::xlab("RV") +
    ggplot2::ylab("Q2 values") +
    theme_graphs
  
  # Return the result
  return(list("Permut_results" = PermRes,
              "Plots" = list("R2val" = p1, 
                             "Q2val" = p2,
                             "R2_and_Q2" = p3,
                             "R2val_hist" = p4, 
                             "Q2val_hist" = p5,
                             "R2_and_Q2_hist" = p6)))
}
