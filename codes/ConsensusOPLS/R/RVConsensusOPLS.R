#' @title RVConsensusOPLS
#' @description
#' Consensus OPLS-DA with RV coefficients weighting and 
#' DQ2 computation for discriminant analysis.
#' 
#' @param data the collection list containing each block of data.
#' @param Y The response.
#' @param maxPcomp Maximum number of Y-predictive components. 
#' @param maxOcomp Maximum number of Y-orthogonal components.
#' @param nfold Number of cross-validation rounds (integer).
#' @param cvType Type of cross-validation used. Either \code{nfold} for n-fold
#' cross-validation, \code{mccv} for Monte Carlo CV or \code{mccvb} for Monte 
#' Carlo class-balanced CV.
#' @param cvFrac numeric. Fraction of observations used in the training set. 
#' Default, 2/3. Fraction of data to be used for \code{mccv} and \code{mccvb}
#' cross-validation
#' @param modelType type of OPLS regression model. Can be defined as \code{reg} 
#' for regression or \code{da} for discriminant analysis. Default \code{da}.
#' @param mc.cores Number of cores for parallel computing. Default: 1.
#' @param verbose Logical which indicates whether the user wants to see the 
#' progress bar printed in the \code{ConsensusOLPSCV} function.
#'
#' @return A consensus OPLS model
#'
#' @examples
#' data(demo_3_Omics)
#' ConsensusOPLS::RVConsensusOPLS(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")], 
#'                                 Y=demo_3_Omics$Y, modelType="da", maxPcomp=1, mc.cores=1, nfold=3)
#' @importFrom parallel mclapply
# #' @keywords internal
#' @export
#' 
RVConsensusOPLS <- function(data,
                            Y,
                            maxPcomp = 1, 
                            maxOcomp = 5, 
                            nfold = 3,
                            cvType = "nfold",
                            cvFrac = 2/3,
                            modelType = "da",
                            mc.cores = 1,
                            verbose = FALSE) {
    # Variable format control
    if (!is.list(data)) stop("data is not a list.")
    if (!is.matrix(Y) && !is.vector(Y) && !is.factor(Y)) stop("Y is not either matrix, vector or factor.")
    if (!is.numeric(maxPcomp)) stop("maxPcomp is not numeric.")
    if (!is.numeric(maxOcomp)) stop("maxOcomp is not numeric.")
    if (!is.numeric(nfold)) stop("nfold is not numeric.")
    if (!is.character(cvType))
        stop("cvType is not a character.")
    else if (!(cvType %in% c("nfold", "mccv", "mccvb")))
        stop("cvType must be `nfold`, `mccv` or `mccvb`.")
    
    if (!is.character(modelType)) stop("modelType is not a character.")
    else if (!(modelType %in% c("reg", "da")))
        stop("modelType must be `reg` or `da`.")
    
    if (!is.logical(verbose)) stop("verbose must be logical `TRUE` or `FALSE`.")
    
    # Check collection dimension
    ntable <- length(data)
    nsample <- nrow(data[[1]])
    nvar <- sapply(data, ncol)
    
    # Initialize parameters
    preProcK <- "mc"
    preProcY <- "mc"
    
    if (modelType == "reg") {
        Y <- as.matrix(Y)
        if (ncol(Y) > 1 || !is.numeric(Y)) stop("modelType is not appropriate to Y.")
        if (all(Y %in% c(0,1))) stop("modelType is preferably `da`.") #TODO: is it possible to do logistic regression?
        koplsScale <- koplsScale(X = Y, centerType = preProcY, scaleType = "no")
        Yc <- koplsScale$X
    } else {
        if (nlevels(as.factor(Y)) > length(Y)/2) { # TODO: something better than this check: if at least 2 samples belong to every class
            stop("modelType is preferably `reg`")
        }
        if (is.vector(Y) || is.factor(Y) || ncol(Y) == 1) {
            Y <- koplsDummy(as.vector(Y))
        }
        if (is.null(colnames(Y))) colnames(Y) <- 1:ncol(Y)
        Yc <- Y
    }
    
    # For each data block
    RA <- mclapply(X = 1:ntable, mc.cores = mc.cores, FUN = function(i) {
        # Produce the kernel of the data block
        tmp <- koplsKernel(X1 = data[[i]], X2 = NULL, Ktype = 'p', params = c(order=1))
        # Frobenius norm of the kernel
        xnorm <- norm(x=tmp, type='F')
        # Normalize the Kernel
        AMat <- tmp/xnorm
        # RV coefficient for AMat
        RV <- (RVmodified(X = AMat, Y = Yc) + 1) / 2
        
        return (list(RV=RV, AMat=AMat))
    })
    names(RA) <- names(data)
    # RV coefficient 
    RV <- mclapply(RA, mc.cores = mc.cores, FUN = function(x) x$RV)
    # Normalized kernel
    AMat <- mclapply(RA, mc.cores = mc.cores, FUN = function(x) x$AMat)
    # calculates the weighted sum of blocks kernel by the RV coeff
    W_mat <- Reduce("+", mclapply(1:ntable, mc.cores = mc.cores, FUN = function(i) {
        RA[[i]]$RV * RA[[i]]$AMat
    }))
    
    # Control maxOcomp
    maxOcomp <- min(c(maxOcomp, nsample, nvar))
    
    # Performs a Kernel-OPLS cross-validation for W_mat
    modelCV <- ConsensusOPLSCV(K = W_mat, Y = Y, maxPcomp = maxPcomp, maxOcomp = maxOcomp, 
                               nfold = nfold, cvType = cvType, preProcK = preProcK, 
                               preProcY = preProcY, cvFrac = cvFrac, 
                               modelType = modelType, verbose = verbose)
    
    Ylarg <- ncol(Y)
    
    # Search for the optimal model based on DQ2
    if (modelType == 'da') {
        mc <- min((maxOcomp+1)*Ylarg, mc.cores)
        mcj <- min(sqrt(mc), Ylarg)
        mci <- max(floor(mc.cores/mcj), 1)
        
        results <- mclapply(X = 0:maxOcomp, mc.cores = mci, FUN = function(i) {
            mclapply(X = 1:Ylarg, mc.cores = mcj, FUN = function(j) {
                # For each Y column, perform the DQ2
                result <- DQ2(Ypred = matrix(data = modelCV$cv$AllYhat[, Ylarg*i+j], #TODO: AllYhat is different from that of Matlab from 3rd column.
                                             ncol = 1), 
                              Y = Y[unlist(modelCV$cv$cvTestIndex), j, drop=F])
                return (result)
            })
        })
        dqq <- do.call(rbind, mclapply(X = 0:maxOcomp, mc.cores = mci, FUN = function(i) {
            do.call(cbind, mclapply(X = 1:Ylarg, mc.cores = mcj, FUN = function(j) {
                return (results[[i+1]][[j]]$dqq) #TODO: why two columns are the same
            }))
        }))
        PRESSD <- do.call(rbind, mclapply(X = 0:maxOcomp, mc.cores = mci, FUN = function(i) {
            do.call(cbind, mclapply(X = 1:Ylarg, mc.cores = mcj, FUN = function(j) {
                return (results[[i+1]][[j]]$PRESSD)
            }))
        }))
        
        dq2 <- rowMeans(dqq)
        index <- maxPcomp + 1 #TODO
        
        # Finds the optimal number of orthogonal components as a function of DQ2
        while (index < (maxOcomp+maxPcomp) && 
               !is.na((dq2[index+1] - dq2[index])) &&
               (dq2[index+1] - dq2[index]) > 0.01) {
            index <- index + 1
        }
        
        # Add DQ2 in the model objects
        modelCV$cv$DQ2Yhat <- dq2 
        # Add optimal number of orthogonal components in the model objects
        modelCV$cv$OrthoLVsOptimalNum <- index - maxPcomp
        
    } else { # if modelType == "reg"
        index <- maxPcomp + 1
        
        # Finds the optimal number of orthogonal components as a function of Q2Yhat
        while (index < (maxOcomp+maxPcomp) && 
               (modelCV$cv$Q2Yhat[index+1] - modelCV$cv$Q2Yhat[index]) > 0.01) {
            index <- index + 1
        }
        # Add optimal number of orthogonal components in the model objects
        modelCV$cv$OrthoLVsOptimalNum <- index - maxPcomp
    }
    
    # Simplifies the name to be used afterwards
    OrthoLVsNum <- modelCV$cv$OrthoLVsOptimalNum
    # print("debug************************")
    # print(OrthoLVsNum)
    # Recompute the optimal model using OrthoLVsNum parameters
    modelCV$Model <- koplsModel(K = W_mat, Y = Y, A = maxPcomp, nox = OrthoLVsNum, 
                                preProcK = preProcK, preProcY = preProcY)
    
    # Adjust Yhat to the selected model size
    #modelCV$cv$Yhat <- modelCV$cv$AllYhat[, Ylarg*OrthoLVsNum + 1:Ylarg, drop=F]

    # Compute the blocks contributions for the selected model
    lambda <- cbind(do.call(rbind,
                            mclapply(X = 1:ntable, mc.cores = mc.cores, 
                                     FUN = function(j) {
                                         diag(crossprod(x = modelCV$Model$scoresP[, 1:maxPcomp], 
                                                        y = crossprod(x = t(AMat[[j]]), 
                                                                      y = modelCV$Model$scoresP[, 1:maxPcomp])))
                                     })),
                    do.call(rbind,
                            mclapply(X = 1:ntable, mc.cores = mc.cores, FUN = function(j) {
                                diag(crossprod(x = modelCV$Model$scoresO[, 1:OrthoLVsNum], 
                                               y = crossprod(x = t(AMat[[j]]), 
                                                             y = modelCV$Model$scoresO[, 1:OrthoLVsNum])))
                            })))
    rownames(lambda) <- names(data)
    colnames(lambda) <- c(paste0("p_", 1:maxPcomp), paste0("o_", 1:OrthoLVsNum))
    
    # Stores raw lambda coefficient values in the model object
    modelCV$Model$lambda_raw <- lambda 
    
    # Normalize the lambda coefficients
    lambda <- sweep(x = lambda, MARGIN = 2, STATS = colSums(lambda), FUN = '/')
    
    # Stores normalized lambda values in the model object
    modelCV$Model$lambda <- lambda 
    
    # Compute the loadings for the selected model size
    loadings.list <- list( 
        mclapply(X = 1:ntable, mc.cores = mc.cores, FUN = function(i) {
            lpi <- tcrossprod(crossprod(x = data[[i]], 
                                        y = modelCV$Model$scoresP[, 1:maxPcomp, drop=F]), 
                              diag(diag(crossprod(modelCV$Model$scoresP[, 1:maxPcomp, drop=F]))^(-1), ncol=maxPcomp))
            colnames(lpi) <- paste0("p", 1:maxPcomp)
            return (lpi)
        }),
        mclapply(X = 1:ntable, mc.cores = mc.cores, FUN = function(i) {
            loi <- tcrossprod(crossprod(x = data[[i]],
                                        y = modelCV$Model$scoresO[, 1:OrthoLVsNum, drop=F]),
                              diag(diag(crossprod(modelCV$Model$scoresO[, 1:OrthoLVsNum, drop=F]))^(-1), ncol=OrthoLVsNum))
            colnames(loi) <- paste0("o", 1:OrthoLVsNum)
            return (loi)
        }))
    loadings <- mclapply(X = 1:ntable, mc.cores = mc.cores, FUN = function(i) { 
        li <- cbind(loadings.list[[1]][[i]], loadings.list[[2]][[i]])
        colnames(li) <- c(paste0("p_", 1:maxPcomp), paste0("o_", 1:OrthoLVsNum))
        return (li)
    })
    names(loadings) <- names(data)
    
    # Add RV coefficients in the model objects
    modelCV$RV <- unlist(RV)
    # Add normalized kernels in the model objects
    modelCV$AMat <- AMat  
    # Add the loadings in the model objects
    modelCV$Model$loadings <- loadings 
    # Add the scores in the model objects
    modelCV$Model$scores <- cbind(modelCV$Model$scoresP, modelCV$Model$scoresO)
    
    # Return the final model
    return (modelCV)
}
