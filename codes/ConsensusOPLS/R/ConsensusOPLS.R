#' @title ConsensusOPLS
#' @description
#' Constructs the consensus OPLS model with an optimal number of orthogonal components
#' for given data blocks and response, and evaluate the model quality w.r.t other models 
#' built with randomly permuted responses.
#' 
#' @param data A list of numeric matrices.
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
#' @param nfold Number of folds performed in n-fold cross-validation. This can be set to the number
#' of samples to perform Leave-One-Out cross validation. Default, 5.
#' @param nMC An integer indicating the number of rounds performed when cvType is `mccv` or `mccvb`.
#' Default, 100.
#' @param cvFrac A numeric value indicating the fraction of observations from \code{data} 
#' used in the training set for `mccv` or `mccvb` cross-validation. Default, 4/5.
#' @param kernelParams List of parameters for the kernel. Default: list(type='p', params = c(order=1.0)).
#' @param mc.cores Number of cores for parallel computing. Default: 1.
#' @param plots A logical indicating if plots are generated. Default: FALSE.
#' @param verbose A logical indicating if the computation progress will be shown. Default, FALSE.
#'
#' @return \code{ConsensusOPLS} returns a list of 
#' \item{\code{optimal}}{ optimal consensus OPLS model}
#' \item{\code{permuted}}{ models with permuted responses}
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
#' @importFrom reshape2 melt
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
                          plots = FALSE,
                          verbose = FALSE) {
    # Variable format control
    if (!is.list(data)) stop("data is not a list.")
    if (!all(sapply(data, is.matrix))) stop("Numeric matrices are required.")
    if (!all(sapply(data, is.numeric))) stop("Numeric matrices are required.")
    if (!is.matrix(Y) && !is.vector(Y) && !is.factor(Y)) stop("Y is not either matrix, vector or factor.")
    if (nperm != as.integer(nperm)) stop("nperm is not an integer.")
    if (maxPcomp != as.integer(maxPcomp)) stop("maxPcomp is not an integer")
    if (maxOcomp != as.integer(maxOcomp)) stop("maxOcomp is not an integer")
    
    if (modelType == "reg") {
        Y <- as.matrix(Y)
        if (ncol(Y) > 1 || !is.numeric(Y)) stop("modelType is not appropriate to Y.")
        if (all(Y %in% c(0,1))) stop("modelType is preferably `da`.") #TODO: is it possible to do logistic regression?
    } else {
        if (nlevels(as.factor(Y)) > length(Y)/2) { # TODO: something better than this check: if at least 2 samples belong to every class
            stop("modelType is preferably `reg`")
        }
        if (is.vector(Y) || is.factor(Y) || ncol(Y) == 1) {
            Y <- koplsDummy(as.vector(Y))
        }
    }
    if (is.null(colnames(Y))) colnames(Y) <- 1:ncol(Y)
    
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
    
    if (plots) {
        # plot blockContribution of the optimal model
        contributions <- perms[[1]]$modelCV$Model$blockContribution
        contributions <- reshape2::melt(contributions)
        colnames(contributions) <- c("block", "comp", "value")
        
        p_contribution <- ggplot(contributions, aes(x=comp, y=value, fill=block)) +
            geom_bar(stat = "identity", position=position_dodge()) +
            theme_light() +
            labs(x = "components", y = "contribution of blocks",
                 fill = "block",
                 title = "contributions")
        
        # plot scores of the optimal model
        scores <- perms[[1]]$modelCV$Model$scores
        scores <- data.frame(scores, response=if (modelType=='da') koplsReDummy(Y) else Y)
        p_scores <- ggplot(scores, aes(x=p_1, y=o_1, colour=response)) +
            geom_point(size=4) +
            labs(x="Predictive",
                 y="Orthogonal",
                 title = "Scores on first orthogonal and predictive components") +
            theme_light()

        # plot loadings of the optimal model
        loadings <- do.call(rbind.data.frame, perms[[1]]$modelCV$Model$loadings)
        loadings$block <- do.call(c, lapply(names(perms[[1]]$modelCV$Model$loadings), function(x) 
            rep(x, nrow(perms[[1]]$modelCV$Model$loadings[[x]]))))
        loadings$variable <- gsub(paste(paste0(names(perms[[1]]$modelCV$Model$loadings), '.'), 
                                        collapse='|'), '', 
                                  rownames(loadings))
        
        p_loadings <- ggplot(loadings, aes(x=p_1, y=o_1, col=block, label=variable)) +
            geom_point(size=2) +
            labs(x="Predictive",
                 y="Orthogonal",
                 title = "Loadings on first orthogonal and predictive components") +
            theme_light()
        
        # plot loadings and VIP of the optimal model
        VIP <- data.frame(VIP = unlist(perms[[1]]$VIP), 
                          variable = unlist(lapply(perms[[1]]$VIP, names)))
        
        loadings_VIP <- merge(loadings, VIP, by="variable")
        loadings_VIP$label <- ifelse(loadings_VIP$VIP > 1, loadings_VIP$variable, NA)
        
        p_vip <- ggplot(loadings_VIP, aes(x=p_1, y=VIP, col=block, label = label)) +
            geom_point(size=2) +
            labs(x="Predictive",
                 y="VIP",
                 title = "VIP versus loadings on predictive components") +
            theme_light()
        
        # plot Q2 permutations
        Q2Yperm <- data.frame(Q2Yperm = permStats$Q2Yhat)
        
        p_q2 <- ggplot(data = Q2Yperm, aes(x = Q2Yperm)) +
            geom_histogram(color="grey", fill="grey") +
            geom_vline(aes(xintercept=permStats$Q2Yhat[1]), color="blue", linetype="dashed", size=1) +
            theme_light() +
            ggtitle("Q2 Permutation test")
        
        # plot DQ2 permutations
        if (modelType=='da') {
            DQ2Yperm <- data.frame(DQ2Yperm = permStats$DQ2Yhat)
            
            p_dq2 <- ggplot(data = DQ2Yperm, aes(x = DQ2Yperm)) +
                geom_histogram(color="grey", fill="grey") +
                geom_vline(aes(xintercept=permStats$DQ2Yhat[1]), color="blue", linetype="dashed", size=1) +
                theme_light() +
                ggtitle("DQ2 Permutation test")
        }
        
        # plot R2Y permutations
        R2Yperm <- data.frame(R2Yperm = permStats$R2Yhat)
        
        p_r2 <- ggplot(data = R2Yperm, aes(x = R2Yperm)) +
            geom_histogram(color="grey", fill="grey") +
            geom_vline(aes(xintercept=permStats$R2Yhat[1]), color="blue", linetype="dashed", size=1) +
            theme_light() +
            ggtitle("R2 Permutation test")
    }
    # Q2Yhat plot
    # hist(x=copls.da$permStats$Q2Yhat, breaks=ceiling(nperm/10), freq = F, 
    #      main="Q2 Permutation test",
    #      xlim=range(copls.da$permStats$Q2Yhat), ylim=c(0,1), xlab='Q2Yhat')
    # lines(density(copls.da$permStats$Q2Yhat), lwd=2)
    # abline(v=copls.da$permStats$Q2Yhat[1], lt='dashed', lwd=2, col='red')
    # text(x=copls.da$permStats$Q2Yhat[1]+diff(range(copls.da$permStats$Q2Yhat))/10, y=0.9,
    #      labels=paste0("p-value=", round(sum(copls.da$permStats$Q2Yhat>=copls.da$permStats$Q2Yhat[1])/(nperm-1), digits=2)))
    # 
    
    # Plotting data
    # permStats_df <- data.frame(
    #     "RV" = c(permStats$RV, 
    #              permStats$RV),
    #     "Values" = c(permStats$R2Yhat, 
    #                  permStats$Q2Yhat),
    #     "Type" = c(rep("R2", times = length(permStats$R2Yhat)),
    #                rep("Q2", times = length(permStats$Q2Yhat)))
    # )
    
    # theme_graphs <- theme_bw() + theme(strip.text = element_text(size=14),
    #                                    axis.title = element_text(size=16),
    #                                    axis.text = element_text(size=14),
    #                                    plot.title = element_text(size=16),
    #                                    legend.title = element_text(size=14))
    # 
    # # Scatterplots
    # p1 <- ggplot2::ggplot(data = permStats_df[which(permStats_df$Type == "RV"), ],
    #                       aes(x = RV, y = Values)) +
    #     ggplot2::geom_point(size = 2.5) +
    #     ggplot2::geom_smooth(method = 'lm') +
    #     ggplot2::xlab("RV") +
    #     ggplot2::ylab("R2 values") +
    #     theme_graphs
    # 
    # p2 <- ggplot2::ggplot(data = permStats_df[which(permStats_df$Type == "Q2"), ],
    #                       aes(x = RV, y = Values)) +
    #     ggplot2::geom_point(size = 2.5) +
    #     ggplot2::geom_smooth(method = 'lm') +
    #     ggplot2::xlab("RV") +
    #     ggplot2::ylab("Q2 values") +
    #     theme_graphs
    # 
    # # Scatterplots all in one
    # p3 <- ggplot2::ggplot(data = permStats_df,
    #                       aes(x = RV, y = Values, col = Type)) +
    #     ggplot2::geom_point(size = 2.5) +
    #     ggplot2::xlab("RV") +
    #     ggplot2::ylab("R2 and Q2 values") +
    #     theme_graphs
    # 
    # # Histograms
    # p4 <- ggplot2::ggplot(data = permStats_df[which(permStats_df$Type == "R2"), ],
    #                       aes(x = Values)) +
    #     ggplot2::geom_histogram(color="black", fill="white") +
    #     ggplot2::geom_density(alpha=.2) +
    #     ggplot2::xlab("Frequency") +
    #     ggplot2::ylab("R2 values") +
    #     theme_graphs
    # 
    # p5 <- ggplot2::ggplot(data = permStats_df[which(permStats_df$Type == "Q2"), ],
    #                       aes(x = Values)) +
    #     ggplot2::geom_histogram(color="black", fill="white") +
    #     ggplot2::geom_density(alpha=.2) +
    #     ggplot2::xlab("Frequency") +
    #     ggplot2::ylab("Q2 values") +
    #     theme_graphs
    # 
    # # Histograms all in one
    # p6 <- ggplot2::ggplot(data = permStats_df,
    #                       aes(x = Values, col = Type)) +
    #     ggplot2::geom_histogram(linewidth = 1) +
    #     ggplot2::geom_density(alpha=.2) +
    #     ggplot2::xlab("RV") +
    #     ggplot2::ylab("Q2 values") +
    #     theme_graphs
    
    return (list(
        optimal = perms[[1]],
        permuted = perms[-1],
        permStats = permStats,
        plots = if (!plots) NULL else list(contribution = p_contribution,
                                           scores       = p_scores,
                                           loadings     = p_loadings,
                                           VIP          = p_vip,
                                           Q2           = p_q2,
                                           DQ2          = if (modelType=='da') p_dq2 else NULL,
                                           R2           = p_r2)))
}
