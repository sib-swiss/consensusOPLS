#' @title \code{ConsensusOPLS} S4 class
#' @description This object presents the result of a Consensus OPLS analysis.
#' @slot nPcomp Number of Y-predictive components (latent variables) of the optimal model.
#' @slot nOcomp Number of Y-orthogonal components (latent variables) of the optimal model.
#' @slot blockContribution Relative contribution of each block (normalized \code{lambda} values) to the
#' latent variables.
#' @slot scores Representation of the samples in the latent variables of the optimal model.
#' @slot loadings Contribution of each block's variables to the latent variables of the optimal
#' model.
#' @slot VIP Variable importance in projection (VIP) for each block of data, assessing the
#' relevance of the variables in explaining the variation in the response.
#' @slot R2X Proportion of variation in data blocks explained by the optimal model.
#' @slot R2Y Proportion of variation in the response explained by the optimal model.
#' @slot Q2 Predictive ability of the optimal model.
#' @slot DQ2 Predictive ability of the optimal model, for discriminant analysis.
#' @slot permStats Q2 and R2Y of models with permuted response.
#' @slot plots Basic plots for the analysis.
#' @slot cv Cross-validation result towards the optimal model.
#' @name ConsensusOPLS-class
#' @rdname ConsensusOPLS-class
#' @import methods utils
# #' @exportMethod generic
#' @export
#' 
setClass("ConsensusOPLS",
         slots = list(
             nPcomp            = "numeric",
             nOcomp            = "numeric",
             blockContribution = "matrix",
             scores            = "matrix",
             loadings          = "list",
             VIP               = "list",
             R2X               = "numeric",
             R2Y               = "numeric",
             Q2                = "numeric",
             DQ2               = "numeric",
             permStats         = "list",
             plots             = "list",
             cv                = "list")
)


#' @noRd
#' 
setMethod("show", "ConsensusOPLS",
          function(object) {
              cat("***Optimal Consensus OPLS model***\n")
              cat("\nNumber of predictive components: ", object@nPcomp, "\n")
              cat("\nNumber of orthogonal components: ", object@nOcomp, "\n")
              cat("\nBlock contribution:\n")
              print(object@blockContribution)
              cat("\nExplained variance R2 in reponse:", object@R2Y[length(object@R2Y)], "\n")
              cat("\nPredictive ability (cross validation Q2): ", object@Q2[length(object@Q2)], "\n")
          }
)


#' @title ConsensusOPLS
#' @description
#' Constructs the consensus OPLS model with an optimal number of orthogonal 
#' components for given data blocks and response, and evaluate the model 
#' quality w.r.t other models built with randomly permuted responses.
#' 
#' @param data A list of data blocks. Each element of the list must be of matrix 
#' type. Rows and columns can have a name, in which case it will be retained 
#' during analysis. Any pre-processing of the data (e.g. scaling) must be 
#' carried out before building the list.
#' @param Y A vector, factor, dummy matrix or numerical matrix for the response. 
#' The type of answer given will condition the model to be used: a numerical 
#' vector for linear regression, a factor or dummy matrix for logistic 
#' regression or a discriminant model.
#' @param maxPcomp Maximum number of Y-predictive components. Default, 1.
#' @param maxOcomp Maximum number of Y-orthogonal components. Default, 5.
#' @param modelType String for type of OPLS regression model, either \code{reg} 
#' for regression or \code{da} for discriminant analysis. Default, \code{da}.
#' @param nperm Number of random permutations desired in response Y. Default, 100.
#' @param cvType String for type of cross-validation used. Either \code{nfold} 
#' for n-fold cross-validation, where \code{nfold} is look up, or \code{mccv} 
#' for Monte Carlo cross-validation, or \code{mccvb} for Monte Carlo 
#' class-balanced cross-validation, where \code{nMC} and \code{cvFrac} are used.
#' Default, \code{nfold}.
#' @param nfold Number of folds performed in n-fold cross-validation. This can 
#' be set to the number of samples to perform Leave-One-Out cross validation. 
#' Default, 5.
#' @param nMC An integer indicating the number of rounds performed when 
#' \code{cvType} is \code{mccv} or \code{mccvb}. Default, 100.
#' @param cvFrac A numeric value indicating the fraction of observations from 
#' \code{data} used in the training set for \code{mccv} or \code{mccvb} 
#' cross-validation. Default, 4/5 = 0.8.
#' @param kernelParams List of parameters for the kernel. Default, 
#' list(type='p', params = c(order=1.0)).
#' @param mc.cores Number of cores for parallel computing. Default, 1.
#' @param plots A logical indicating if plots are generated. For more aesthetic 
#' graphics, please refer to the package thumbnail. Interpretation help is also 
#' provided. Default, FALSE.
#' @param verbose A logical indicating if detailed information (cross
#' validation) will be shown. Default, FALSE.
#'
#' @return An object of class \code{ConsensusOPLS}.
#' @examples
#' data(demo_3_Omics)
#' res <- ConsensusOPLS(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")], 
#'                      Y=demo_3_Omics$Y,
#'                      maxPcomp=1, maxOcomp=2, 
#'                      modelType='da',
#'                      nperm=5)
#' str(res)
#' @importFrom reshape2 melt
#' @import utils ggplot2 parallel 
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
    if (!is.matrix(Y) && !is.vector(Y) && !is.factor(Y)) stop("Y is not either matrix, vector or factor.")
    if (nperm != as.integer(nperm)) stop("nperm is not an integer.")
    if (maxPcomp != as.integer(maxPcomp)) stop("maxPcomp is not an integer")
    if (maxOcomp != as.integer(maxOcomp)) stop("maxOcomp is not an integer")
    
    if (modelType == "reg") {
        if (is.factor(Y)) {
            warning("Logistic regression is performed.")
            if (nlevels(Y) > 2) stop("Y should have two levels.")
            Y <- koplsDummy(as.vector(Y))[, 1, drop=F]
        } else {
            Y <- as.matrix(Y)
        }
        if (ncol(Y) > 1 || !is.numeric(Y)) stop("modelType is preferably `da`.")
    } else {
        if (nlevels(as.factor(Y)) > length(Y)/2) { # TODO: something better than this check: if at least 2 samples belong to every class
            stop("modelType is preferably `reg`.")
        }
        if (is.vector(Y) || is.factor(Y) || ncol(Y) == 1) {
            Y <- koplsDummy(as.vector(Y))
        }
    }
    
    if (is.matrix(Y) && is.null(colnames(Y))) colnames(Y) <- 1:ncol(Y)
    
    # Init parameters
    permStats <- list()
    # Loaded packages
    loaded.package.names <- c(
        sessionInfo()$basePkgs,
        names(sessionInfo()$otherPkgs))
    
    # Create parallel cluster
    globLibPaths <- .libPaths()
    cl <- makeCluster(mc.cores)
    clusterExport(cl,
                  ls(all.names=TRUE, envir=globalenv()),
                  envir=.GlobalEnv)
    # Load the packages on all the cluster
    parLapply(cl, 1:length(cl), function(i) {
        lapply(loaded.package.names, function(lpn) {
            .libPaths(globLibPaths)
            library(lpn, character.only=TRUE)
        })
    })
    
    # Permutations
    #perms <- mclapply(X=1:(1+nperm), mc.cores=mc.cores, function(i) {
    perms <- parLapply(cl, X=1:(1+nperm), function(i) {
        # Fix the random seed
        set.seed(i)
        
        # Random permutation of Y rows
        if (i==1) {
            Ys <- Y
        } else {
            Ys <- Y[sample(x = 1:nrow(Y), size = nrow(Y), replace = FALSE, prob = NULL), , drop=F]
        }
            
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
    #permStats$lvnum  <- unlist(mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
    permStats$lvnum  <- unlist(parLapply(cl, X=1:(1+nperm), function(i) {
        perms[[i]]$modelCV$cv$nOcompOpt + maxPcomp
    }))
    #permStats$R2Yhat  <- unlist(mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
    permStats$R2Yhat  <- unlist(parLapply(cl, X=1:(1+nperm), function(i) {
        utils::tail(perms[[i]]$modelCV$Model$R2Yhat, 1)
    }))
    #permStats$DQ2Yhat <- unlist(mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
    permStats$DQ2Yhat <- unlist(parLapply(cl, X=1:(1+nperm), function(i) {
        perms[[i]]$modelCV$cv$DQ2Yhat[permStats$lvnum[i]]
    }))
    #permStats$Q2Yhat  <- unlist(mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
    permStats$Q2Yhat  <- unlist(parLapply(cl, X=1:(1+nperm), function(i) {
        perms[[i]]$modelCV$cv$Q2Yhat[permStats$lvnum[i]]
    }))
    # permStats$PredAc <- unlist(mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
    #     perms[[i]]$modelCV$da$tot_sens[2]
    # }))
    #permStats$Y      <- mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
    permStats$Y      <- unlist(parLapply(cl, X=1:(1+nperm), function(i) {
        perms[[i]]$Ys
    }))
    #permStats$RV     <- unlist(mclapply(1:(1+nperm), mc.cores=mc.cores, function(i) {
    permStats$RV     <- unlist(parLapply(cl, X=1:(1+nperm), function(i) {
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
            labs(x = "Predictive",
                 y = "Orthogonal",
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
            labs(x = "Predictive",
                 y = "Orthogonal",
                 title = "Loadings on first orthogonal and predictive components") +
            theme_light()
        
        # plot loadings and VIP of the optimal model
        VIP <- data.frame(VIP = unlist(perms[[1]]$VIP), 
                          variable = unlist(lapply(perms[[1]]$VIP, names)))
        
        loadings_VIP <- merge(loadings, VIP, by="variable")
        loadings_VIP$label <- ifelse(loadings_VIP$VIP > 1, loadings_VIP$variable, NA)
        
        p_vip <- ggplot(loadings_VIP, aes(x=p_1, y=VIP, col=block, label = label)) +
            geom_point(size=2) +
            labs(x = "Predictive",
                 y = "VIP",
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
    
    ## Stop parallel clusters
    stopCluster(cl)
    
    return (new("ConsensusOPLS",
                nPcomp            = perms[[1]]$modelCV$Model$params$nPcomp,
                nOcomp            = perms[[1]]$modelCV$Model$params$nOcomp,
                blockContribution = perms[[1]]$modelCV$Model$blockContribution,
                scores            = perms[[1]]$modelCV$Model$scores,
                loadings          = perms[[1]]$modelCV$Model$loadings,
                VIP               = perms[[1]]$VIP,
                R2X               = perms[[1]]$modelCV$Model$R2X,
                R2Y               = perms[[1]]$modelCV$Model$R2Yhat,
                Q2                = perms[[1]]$modelCV$cv$Q2Yhat[1:(perms[[1]]$modelCV$Model$params$nOcomp+1)],
                DQ2               = if (modelType=='da') 
                    perms[[1]]$modelCV$cv$DQ2Yhat[1:(perms[[1]]$modelCV$Model$params$nOcomp+1)] else numeric(),
                permStats         = list(Q2Y=permStats$Q2Yhat,
                                         DQ2Y=permStats$DQ2Yhat,
                                         R2Y=permStats$R2Yhat),
                plots             = if (!plots) list() else 
                    list(contribution = p_contribution,
                         scores       = p_scores,
                         loadings     = p_loadings,
                         VIP          = p_vip,
                         Q2           = p_q2,
                         DQ2          = if (modelType=='da') p_dq2 else NULL,
                         R2           = p_r2),
                cv                = if (verbose) perms[[1]]$modelCV$cv else list())
            )
}
