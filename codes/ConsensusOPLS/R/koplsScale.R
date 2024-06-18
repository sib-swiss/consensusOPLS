#' @title koplsScale
#' @description Function for mean-centering and scaling of a matrix.
#'
#' @param X matrix. The matrix to be centered and scaled.
#' @param centerType character. Indicates the centering method of the \code{X} 
#' matrix: mean centering (\code{mc}) or no centering (\code{no}). Default is 
#' \code{no} centering.
#' @param scaleType character. Indicates the scaling method of the \code{X} 
#' matrix: \code{uv} for unit variance scaling (\code{pa}) for Pareto scaling or 
#' \code{no} for no scaling. Default is \code{no} scaling.
#'
#' @returns A list containing the following entries:
#' \item{centerType}{ character. Indicates the centering method of the X matrix.}
#' \item{scaleType}{ character. Indicates the scaling method of the X matrix.}
#' \item{meanV}{ vector. Contains the mean values for all columns in X.}
#' \item{stdV}{ vector. Contains the standard deviations for all columns in X.}
#' \item{matrix}{ matrix. Original input matrix X, scaled according to 
#' \code{centerType} and \code{scaleType}.}
#' 
#' @examples
#' X <- matrix(data = c(1,4,7, 8,4,0, 3,6,9), nrow = 3)
#' Y <- ConsensusOPLS:::koplsScale(X = X, centerType = "mc", scaleType = "pa")
#' Y$X
#' 
#' @import stats
#' @keywords internal
#' @noRd
#' 
koplsScale <- function(X, centerType = "no", scaleType = "no"){
    # Variable format control
    if (!is.matrix(X)) stop("X is not a matrix.")
    if (!is.character(centerType)) stop("centerType is not a character.")
    else if(!(centerType %in% c("mc", "no")))
        stop("centerType must be `mc` or `no`.")
    
    if (!is.character(scaleType)) stop("scaleType is not a character.")
    else if(!(scaleType %in% c("uv", "pa", "no")))
        stop("scaleType must be `uv`, `pa` or `no`.")
    
    # Calculation of dispersion parameters before center and scale matrix
    meanV <- colMeans(X)
    stdV <- apply(X = X, MARGIN = 2, FUN = function(X) stats::sd(X))
    
    # Center the matrix
    if (centerType == "mc")
        X <- scale(x = X, center = TRUE, scale = FALSE)
    
    # Scale the matrix
    if (scaleType == "uv")
        X <- scale(x = X, center = FALSE, scale = TRUE)
    else if (scaleType == "pa")
        X <- apply(X = X, MARGIN = 2, FUN = function(col){ col/sqrt(stats::sd(col))})
    
    # Return a list with all parameters
    return(list("centerType" = centerType,
                "scaleType" = scaleType,
                "meanV" = meanV,
                "stdV" = stdV,
                "X" = X))
}



#' @title koplsRescale
#' @description Scales a matrix based on pre-defined parameters from a scaling
#' object defined in a list (result of \code{koplsScale} function).
#' 
#' @param scaleS list. It contains scaling parameters.
#' @param varargin matrix. If defined, this matrix will be scaled and returned.
#' Otherwise the original data set in the \code{scaleS} object will be scaled 
#' and returned. 
#'
#' @returns A list containing the following entries:
#' \item{centerType}{ character. Indicates the centering method of the X matrix:
#' \code{mc} (mean-centering) or \code{no} (no centering).}
#' \item{scaleType}{ character. Indicates the scaling method of the X matrix: 
#' \code{uv} (unit variance), \code{pa} (pareto) or \code{no} (no scaling).}
#' \item{meanV}{ vector. Contains the mean values for all columns in X.}
#' \item{stdV}{ vector. Contains the standard deviations for all columns in X.}
#' \item{X}{ matrix. Scaled version of \code{varargin}, if defined. 
#' Otherwise, scaled version of \code{scaleS$X} from input. Scaling is done 
#' according to \code{centerType} and \code{scaleType}.
#' }
#' 
#' @examples
#' data <- matrix(data = c(-1.732051, 0, 1.732051, 
#'                         2, 0,-2,
#'                         -1.732051, 0, 1.732051),
#'                nrow = 3, ncol = 3)
#' scaleS <- list("centerType" = "mc", "scaleType" = "pa", "meanV" = 0, 
#'                "stdV" = 1.581139, "matrix" = data)
#' test <- ConsensusOPLS:::koplsRescale(scaleS = scaleS, varargin = NULL)
#' test
#' test$X
#' 
#' @keywords internal
#' @noRd
#' 
koplsRescale <- function(scaleS, varargin = NULL){
    # Variable format control
    if (!is.list(scaleS)) stop("scaleS must be a list (result of `koplsScale()`).")
    if (!is.null(varargin)) {
        if (!is.matrix(varargin)) stop("varargin must be a matrix.")
        X <- varargin
    } else {
        X <- scaleS$X
    }
    
    # Center the matrix
    if (scaleS$centerType == "mc") 
        X <- X + scaleS$meanV
    
    # Scale the matrix
    if (scaleS$scaleType == "uv") 
        X <- X * scaleS$stdV
    if (scaleS$scaleType == "pa") 
        X <- X * sqrt(scaleS$stdV)
    
    return (list("centerType" = "no",
                 "scaleType" = "no",
                 "meanV" = scaleS$meanV, 
                 "stdV" = scaleS$stdV,
                 "X" = X))
}



#' @title koplsScaleApply
#' @description Applies scaling from external scaling objects on a matrix X.
#' 
#' @param model list. An object containing scaling parameters.
#' @param X matrix. The matrix to be scaled according to model parameters.
#'
#' @returns A list containing the following entries:
#' \item{centerType}{ character. Indicates the centering method of the X matrix:
#' \code{mc} (mean-centering) or \code{no} (no centering).}
#' \item{scaleType}{ character. Indicates the scaling method of the X matrix: 
#' \code{uv} (unit variance), \code{pa} (pareto) or \code{no} (no scaling).}
#' \item{meanV}{ vector. Contains the mean values for all columns in X.}
#' \item{stdV}{ vector. Contains the standard deviations for all columns in X.}
#' \item{X}{ matrix. Scaled version of \code{varargin}, if defined. 
#' Otherwise, scaled version of \code{scaleS$X} from input. Scaling is done 
#' according to \code{centerType} and \code{scaleType}.
#' }
#' 
#' @examples
#' X <- matrix(data = c(1,4,7, 8,4,0, 3,6,9), nrow = 3)
#' Y <- ConsensusOPLS:::koplsScale(X = X, centerType = "mc", scaleType = "pa")
#' Z <- ConsensusOPLS:::koplsScaleApply(model = Y, X = Y$X)
#' Z$X
#' 
#' @keywords internal
#' @noRd
#' 
koplsScaleApply <- function(model, X){
    # Variable format control
    if (!is.list(model)) stop("model is not a list with scaling parameters.")
    if (!is.matrix(X)) stop("X is not a matrix.")
    
    # Center the matrix
    if (model$centerType == "mc") X <- X - model$meanV
    
    # Scale the matrix
    if (model$scaleType == "uv") X <- X / model$stdV
    if (model$scaleType == "pa") X <- X / sqrt(model$stdV)
    
    # Return a list with all parameters
    return (list("centerType" = model$centerType,
                 "scaleType" = model$scaleType,
                 "meanV" = model$meanV,
                 "stdV" = model$stdV,
                 "X" = X))
}
