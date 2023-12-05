#' @title koplsSensSpec
#' @description Calculates sensitivity and specificity in a class-wise fashion.
#'
#' @param trueClass matrix. Row vector of true class assignments (template). 
#' @param predClass matrix. Matrix (or row vector) of class assignments to be 
#' compared.
#'
#' @return 
#' \item{data}{ A data frame containing, for each class (on lines), the number 
#' of True Positives \code{TP}, True Negatives \code{TN}, False Positives 
#' \code{FP} and False Negatives \code{FN}. It also contains (still in columns), 
#' the number of subjects \code{N}, the sensitivity \code{sens} and the 
#' specificity \code{spec}. The last two columns correspond to the mean 
#' sensitivity \code{meanSens} and mean specificity \code{meanSpec} (which are 
#' equals to the value for each class). The last line shows the total values for 
#' each column (named \code{tot}).}
#' 
#' @examples
#' trueClass <- sample(x = 2:6, size = 100, replace = TRUE)
#' predClass <- sample(x = 1:5, size = 100, replace = TRUE)
#' test <- ConsensusOPLS:::koplsSensSpec(trueClass = trueClass, 
#'                                       predClass = predClass)
#' test
#' 
#' @keywords internal

koplsSensSpec <- function(trueClass, predClass) {
    # Variable format control
    if (!is.data.frame(trueClass) & !is.matrix(trueClass)) {
        warning("trueClass is neither a matrix nor a data.frame,
            so it was converted into a matrix.")
        trueClass <- as.matrix(trueClass)
    }
    if (!is.data.frame(predClass) & !is.matrix(predClass)) {
        warning("predClass is neither a matrix nor a data.frame, 
            so it was converted into a matrix.")
        predClass <- as.matrix(predClass)
    }
    
    # Check dimensions
    if (is.vector(trueClass) & is.vector(predClass)) {
        if (length(predClass) != length(trueClass)) {
            stop("Template vector (trueClass) differs in length from the vector 
           (predClass) to be compared.")
        }
        
        # Forces matrix conversion for the code below
        predClass <- as.matrix(predClass)
    } else {
        if (is.vector(trueClass) & is.matrix(predClass)) {
            if (nrow(predClass) != length(trueClass)) {
                stop("Template vector (trueClass) differs in length from the matrix 
             (predClass) to be compared.")
            }
        }
    }
    
    # Contingency table
    contingency <- table(trueClass, predClass)
    label_class <- as.character(sort(union(x = trueClass, y = predClass)))
    
    # Define TruePositifs
    TP <- sapply(label_class,
                 FUN = function(x) {
                     if (x %in% rownames(contingency) && x %in% colnames(contingency)) {
                         contingency[x, x]
                     } else{0}
                 })
    # Define TrueNegatifs
    TN <- sapply(label_class,
                 FUN = function(x) {
                     sum(TP) - if(x %in% rownames(contingency) && x %in% colnames(contingency)){
                         contingency[x, x]
                     } else{0}
                 })
    # Define FalsePositifs
    FN <- sapply(label_class,
                 FUN = function(x) {
                     if(x %in% rownames(contingency)) {
                         sum(contingency[x, colnames(contingency) != x])
                     } else 0
                 })
    # Define FalseNegatifs
    FP <- sapply(label_class,
                 FUN = function(x) {
                     if (x %in% colnames(contingency)) {
                         sum(contingency[rownames(contingency) != x, x])
                     } else 0
                 })
    
    # Sensitivity for each class
    sens <- TP / (TP + FN)
    # Specificity for each class
    spec <- TN / (TN + FP)
    
    # Correct missing values
    sens[is.na(sens)] <- 0
    spec[is.na(spec)] <- 0
    
    # Results for each class
    results <- data.frame(TP = c(TP, sum(TP)),
                          TN = c(TN, sum(TN)),
                          FP = c(FP, sum(FP)),
                          FN = c(FN, sum(FN)),
                          N = c(TP+TN+FP+FN, sum(TP+FP)),
                          sens = c(sens, sum(TP)/ sum(TP + FN)),
                          spec = c(spec, sum(TN)/ sum(TN + FP)),
                          meanSens = c(sens, mean(sens)),
                          meanSpec = c(spec, mean(spec)))
    
    # Change rownames
    rownames(results) <- c(label_class, "tot")
    
    #Return the data frame
    return (results)
}
