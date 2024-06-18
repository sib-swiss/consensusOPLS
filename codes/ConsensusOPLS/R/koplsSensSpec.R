#' @title koplsSensSpec
#' @description Calculates sensitivity and specificity in a class-wise fashion.
#'
#' @param trueClass matrix. Row vector of true class assignments (template). 
#' @param predClass matrix. Matrix (or row vector) of class assignments to be 
#' compared.
#' @param labelClass A vector of all potential classes. Default, NULL, labels found
#' in \code{trueClass} and \code{codeClass} are considered.
#'
#' @returns
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
#' @noRd
#' 
koplsSensSpec <- function(trueClass, predClass, labelClass = NULL) {
    # Contingency table
    contingency <- table(trueClass, as.character(predClass))
    if (is.null(labelClass))
        labelClass <- as.character(sort(x = union(x = trueClass, y = predClass)))
    
    # Define TruePositifs
    TP <- sapply(labelClass,
                 FUN = function(x) {
                     if (x %in% rownames(contingency) && 
                         x %in% colnames(contingency)) {
                         contingency[x, x]
                     } else 0
                 })
    # Define TrueNegatifs
    TN <- sapply(labelClass,
                 FUN = function(x) {
                     sum(TP) - if(x %in% rownames(contingency) && 
                                  x %in% colnames(contingency)){
                         contingency[x, x]
                     } else 0
                 })
    # Define FalsePositifs
    FN <- sapply(labelClass,
                 FUN = function(x) {
                     if(x %in% rownames(contingency)) {
                         sum(contingency[x, colnames(contingency) != x])
                     } else 0
                 })
    # Define FalseNegatifs
    FP <- sapply(labelClass,
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
    # results <- data.frame(TP = c(TP, sum(TP)),
    #                       TN = c(TN, sum(TN)),
    #                       FP = c(FP, sum(FP)),
    #                       FN = c(FN, sum(FN)),
    #                       N = c(TP+TN+FP+FN, sum(TP+FP)),
    #                       sens = c(sens, sum(TP)/ sum(TP + FN)),
    #                       spec = c(spec, sum(TN)/ sum(TN + FP)),
    #                       meanSens = c(sens, mean(sens)),
    #                       meanSpec = c(spec, mean(spec)))
    results <- data.frame(TP=TP,
                          TN=TN,
                          FP=FP,
                          FN=FN,
                          sens=sens,
                          spec=spec)
    
    return (results)
}
