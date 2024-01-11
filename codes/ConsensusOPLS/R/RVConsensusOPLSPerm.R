#' @title RVConsensusOPLSPerm
#' @description
#' Y random permutations for model validity estimation, with leave-on-out CV.
#' 
#' @param data list. The collection list containing each block of data.
#' @param Y matrix. The response matrix.
#' @param nbruns numeric. Number of random permutations. 
#' @param PredLVs numeric. Number of predictive components.
#' @param maxOrtholvs numeric. Maximum number of orthogonal LVs to compute.
#' @param modelType type of OPLS regression model. Can be defined as \code{reg} 
#' for regression or \code{da} for discriminant analysis. Default value is
#' \code{da}.
#' @param cvType Type of cross-validation used. Either \code{nfold} for n-fold
#' @param mc.cores Number of cores for parallel computing. Default: 1.
#'
#' @return 
#' A list with the two results plots: \code{plot_R2val} for RV vs R2val, \code{plot_Q2val} for RV vs Q2val.
#'
#' @examples
#' data(demo_3_Omics)
#' RVConsensusOPLSPerm(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")], 
#'                     Y=demo_3_Omics$Y, nbruns=5, PredLVs=1, maxOrtholvs=2, modelType = 'da')
#' @importFrom utils tail
#' @importFrom parallel mclapply
#' @import ggplot2
#' @export
#' 
RVConsensusOPLSPerm <- function(data, Y, nbruns, PredLVs, maxOrtholvs, modelType = 'da', cvType = 'nfold', mc.cores = 1) {
    # Variable format control
    if (!is.list(data)) stop("data is not a list.")
    if (!is.matrix(Y)) stop("Y is not a matrix.")
    if (!is.numeric(nbruns)) stop("A is not numeric.")
    if (!is.numeric(PredLVs)) stop("PredLVs is not numeric.")
    if (!is.numeric(maxOrtholvs)) stop("maxOrtholvs is not numeric.")
    
    # Init parameters
    NbObs <- nrow(data[[1]])
    PermRes <- list()
    
    # Permutations
    perms <- mclapply(X=1:(1+nbruns), mc.cores=mc.cores, function(i) {
        # Fix the random seed
        set.seed(i)
        
        # Random permutation of Y rows
        if (i==1) 
            Ys <- Y
        else 
            Ys <- Y[sample(x = 1:nrow(Y), size = nrow(Y), replace = FALSE, prob = NULL), ]
        
        # Redo the Consensus OPLS-DA with RV coefficients weighting
        modelCV <- RVConsensusOPLS(data = data, Y = Ys, A = PredLVs, 
                                   maxOrtholvs = maxOrtholvs, nrcv = NbObs,
                                   cvType = cvType, modelType = modelType, 
                                   mc.cores = 1,
                                   verbose = FALSE)
        return (list(Ys=Ys,
                     modelCV=modelCV))
    })
    PermRes$lvnum <- unlist(mclapply(1:(1+nbruns), mc.cores=mc.cores, function(i) {
        perms[[i]]$modelCV$cv$OrthoLVsOptimalNum + PredLVs
    }))
    PermRes$R2val <- unlist(mclapply(1:(1+nbruns), mc.cores=mc.cores, function(i) {
        tail(perms[[i]]$modelCV$Model$R2Yhat, 1)
    }))
    PermRes$DQ2val <- unlist(mclapply(1:(1+nbruns), mc.cores=mc.cores, function(i) {
        perms[[i]]$modelCV$cv$DQ2Yhat[PermRes$lvnum[i]]
    }))
    PermRes$Q2val <- unlist(mclapply(1:(1+nbruns), mc.cores=mc.cores, function(i) {
        perms[[i]]$modelCV$cv$Q2Yhat[PermRes$lvnum[i]]
    }))
    PermRes$PredAc <- unlist(mclapply(1:(1+nbruns), mc.cores=mc.cores, function(i) {
        perms[[i]]$modelCV$da$tot_sens[2]
    }))
    PermRes$Y <- mclapply(1:(1+nbruns), mc.cores=mc.cores, function(i) {
        perms[[i]]$Ys
    })
    PermRes$RV <- unlist(mclapply(1:(1+nbruns), mc.cores=mc.cores, function(i) {
        RVmodified(X = Y, Y = perms[[i]]$Ys)
    }))
    
    # Plotting data
    PermRes_df <- data.frame(
        "RV" = c(PermRes$RV, 
                 PermRes$RV),
        "Values" = c(PermRes$R2val, 
                     PermRes$Q2val),
        "Type" = c(rep("RV", times = length(PermRes$R2val)),
                   rep("Q2", times = length(PermRes$Q2val)))
    )
    
    theme_graphs <- theme_bw() + theme(strip.text = element_text(size=14),
                                       axis.title = element_text(size=16),
                                       axis.text = element_text(size=14),
                                       plot.title = element_text(size=16),
                                       legend.title = element_text(size=14))
    
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
        ggplot2::geom_histogram(linewidth = 1)+
        ggplot2::geom_density(alpha=.2)+
        ggplot2::xlab("RV") +
        ggplot2::ylab("Q2 values") +
        theme_graphs
    
    return (list("PermRes" = PermRes,
                 "Plots" = list("R2val" = p1, 
                                "Q2val" = p2,
                                "R2_and_Q2" = p3,
                                "R2val_hist" = p4, 
                                "Q2val_hist" = p5,
                                "R2_and_Q2_hist" = p6)))
}
