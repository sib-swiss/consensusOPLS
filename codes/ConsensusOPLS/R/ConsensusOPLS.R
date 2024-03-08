#' @title ConsensusOPLS
#' @description
#' Constructs the consensus OPLS model with an optimal number of orthogonal components
#' for given data blocks and response, and evaluate the model quality w.r.t other models 
#' built with randomly permuted responses.
#' 
#' @param data A list of data blocks.
#' @param Y A vector, factor, dummy matrix or numerical matrix for the response.
#' @param maxPcomp Maximum number of Y-predictive components. Default, 1.
#' @param maxOcomp Maximum number of Y-orthogonal components. Default, 5.
#' @param modelType Type of OPLS regression model, either `reg` for regression or 
#' `da` for discriminant analysis. Default, `da`.
#' @param nperm Number of random permutations. Default, 100.
#' @param cvType Type of cross-validation used. Either `nfold` for n-fold cross-validation, 
#' where \code{nfold} is look up, or `mccv` for Monte Carlo CV, or `mccvb` 
#' for Monte Carlo class-balanced cross-validation, where \code{nMC} and \code{cvFrac} are used.
#' Default, `nfold`.
#' @param nfold Number of folds performed in n-fold cross-validation. Default, 5.
#' @param nMC An integer indicating the number of rounds performed when cvType is `mccv` or `mccvb`.
#' Default, 100.
#' @param cvFrac A numeric value indicating the fraction of observations from \code{data} 
#' used in the training set for `mccv` or `mccvb` cross-validation. Default, 4/5.
#' @param mc.cores Number of cores for parallel computing. Default: 1.
#' @param kernelParams List of parameters for the kernel. Default: list(type='p', params = c(order=1.0)).
#' @param verbose A logical indicating if the computation progress will be shown. Default, FALSE.
#'
#' @return \code{ConsensusOPLS} returns a list of 
#' \item{\code{optModel}}{ optimal consensus OPLS model}
#' \item{\code{permModels}}{ models with permuted responses}
#' \item{\code{permStats}}{ permutation statistics}
#' \item{\code{plots}}{ plots}
#'
#' @examples
#' data(demo_3_Omics)
#' res <- ConsensusOPLS(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")], 
#'                      Y=demo_3_Omics$Y,
#'                      maxPcomp=1, maxOcomp=2, 
#'                      modelType='da',
#'                      nperm=5)
#' str(res)
#' @importFrom parallel mclapply
#' @import ggplot2
#' @export
#' 
ConsensusOPLS <- function(data,
                          Y,
                          maxPcomp = 1,
                          maxOcomp = 5,
                          modelType = 'da',
                          nperm = 100,
                          cvType = 'nfold',
                          nfold = 5,
                          nMC = 100,
                          cvFrac = 4/5,
                          kernelParams = list(type='p', params = c(order=1.0)),
                          mc.cores = 1,
                          verbose = FALSE) {
    # Variable format control
    if (!is.list(data)) stop("data is not a list.")
    if (!is.matrix(Y) && !is.vector(Y) && !is.factor(Y)) stop("Y is not either matrix, vector or factor.")
    if (nperm != as.integer(nperm)) stop("nperm is not an integer.")
    if (maxPcomp != as.integer(maxPcomp)) stop("maxPcomp is not an integer")
    if (maxOcomp != as.integer(maxOcomp)) stop("maxOcomp is not an integer")
    
    # Init parameters
    permStats <- list()
    
    # Permutations
    perms <- mclapply(X=1:(1+nperm), mc.cores=mc.cores, function(i) {
        # Fix the random seed
        set.seed(i)
        
        # Random permutation of Y rows
        if (i==1) 
            Ys <- Y
        else 
            Ys <- Y[sample(x = 1:nrow(Y), size = nrow(Y), replace = FALSE, prob = NULL), , drop=F]
        
        # Redo the Consensus OPLS-DA with RV coefficients weighting
        modelCV <- RVConsensusOPLS(data = data,
                                   Y = Ys,
                                   maxPcomp = maxPcomp,
                                   maxOcomp = maxOcomp,
                                   modelType = modelType,
                                   cvType = cvType,
                                   nfold = nfold,
                                   nMC = nMC,
                                   cvFrac = cvFrac,
                                   mc.cores = 1,
                                   kernelParams = kernelParams,
                                   verbose = verbose)
        VIP <- VIP(data = data, Y = Ys, model=modelCV$Model)
        
        return (list(Ys=Ys,
                     modelCV=modelCV,
                     VIP=VIP)
        )
    })
    permStats$lvnum  <- unlist(mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
        perms[[i]]$modelCV$cv$nOcompOpt + maxPcomp
    }))
    permStats$R2Yhat  <- unlist(mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
        utils::tail(perms[[i]]$modelCV$Model$R2Yhat, 1)
    }))
    permStats$DQ2Yhat <- unlist(mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
        perms[[i]]$modelCV$cv$DQ2Yhat[permStats$lvnum[i]]
    }))
    permStats$Q2Yhat  <- unlist(mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
        perms[[i]]$modelCV$cv$Q2Yhat[permStats$lvnum[i]]
    }))
    # permStats$PredAc <- unlist(mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
    #     perms[[i]]$modelCV$da$tot_sens[2]
    # }))
    permStats$Y      <- mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
        perms[[i]]$Ys
    })
    permStats$RV     <- unlist(mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
        RVmodified(X = Y, Y = perms[[i]]$Ys)
    }))
    
    # Plotting data
    permStats_df <- data.frame(
        "RV" = c(permStats$RV, 
                 permStats$RV),
        "Values" = c(permStats$R2Yhat, 
                     permStats$Q2Yhat),
        "Type" = c(rep("R2", times = length(permStats$R2Yhat)),
                   rep("Q2", times = length(permStats$Q2Yhat)))
    )
    
    theme_graphs <- theme_bw() + theme(strip.text = element_text(size=14),
                                       axis.title = element_text(size=16),
                                       axis.text = element_text(size=14),
                                       plot.title = element_text(size=16),
                                       legend.title = element_text(size=14))
    
    # Scatterplots
    p1 <- ggplot2::ggplot(data = permStats_df[which(permStats_df$Type == "RV"), ],
                          aes(x = RV, y = Values)) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::geom_smooth(method = 'lm') +
        ggplot2::xlab("RV") +
        ggplot2::ylab("R2 values") +
        theme_graphs
    
    p2 <- ggplot2::ggplot(data = permStats_df[which(permStats_df$Type == "Q2"), ],
                          aes(x = RV, y = Values)) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::geom_smooth(method = 'lm') +
        ggplot2::xlab("RV") +
        ggplot2::ylab("Q2 values") +
        theme_graphs
    
    # Scatterplots all in one
    p3 <- ggplot2::ggplot(data = permStats_df,
                          aes(x = RV, y = Values, col = Type)) +
        ggplot2::geom_point(size = 2.5) +
        ggplot2::xlab("RV") +
        ggplot2::ylab("R2 and Q2 values") +
        theme_graphs
    
    # Histograms
    p4 <- ggplot2::ggplot(data = permStats_df[which(permStats_df$Type == "R2"), ],
                          aes(x = Values)) +
        ggplot2::geom_histogram(color="black", fill="white") +
        ggplot2::geom_density(alpha=.2) +
        ggplot2::xlab("Frequency") +
        ggplot2::ylab("R2 values") +
        theme_graphs
    
    p5 <- ggplot2::ggplot(data = permStats_df[which(permStats_df$Type == "Q2"), ],
                          aes(x = Values)) +
        ggplot2::geom_histogram(color="black", fill="white") +
        ggplot2::geom_density(alpha=.2) +
        ggplot2::xlab("Frequency") +
        ggplot2::ylab("Q2 values") +
        theme_graphs
    
    # Histograms all in one
    p6 <- ggplot2::ggplot(data = permStats_df,
                          aes(x = Values, col = Type)) +
        ggplot2::geom_histogram(linewidth = 1) +
        ggplot2::geom_density(alpha=.2) +
        ggplot2::xlab("RV") +
        ggplot2::ylab("Q2 values") +
        theme_graphs
    
    return (list(
        optimal = perms[[1]],
        permuted = perms[-1],
        permStats = permStats,
        plots = list("R2Yhat" = p1, 
                     "Q2Yhat" = p2,
                     "R2_and_Q2" = p3,
                     "R2Yhat_hist" = p4, 
                     "Q2Yhat_hist" = p5,
                     "R2_and_Q2_hist" = p6)))
}
