#' @title DQ2
#' @description
#' Calculate the discriminant Q2 value as an improvement for the Q2 value. DQ2
#' does not penalize class predictions beyond the class label value, i.e. 
#' predicted values < 0 for 0-class and predicted values > 1 for 1-class.
#' Westerhuis, J.A., van Velzen, E.J.J., Hoefsloot, H.C.J. et al. Discriminant
#' Q2 (DQ2) for improved discrimination in PLSDA models. Metabolomics 4,
#' 293–296 (2008). https://doi.org/10.1007/s11306-008-0126-2
#'
#' @param Ypred A vector of predicted values.
#' @param Y A vector of two-class 0/1 labels.
#'
#' @returns
#' \item{DQ2}{numeric. Discriminant Q2 value.}
#' \item{PRESSD}{matrix. The adjusted prediction error sum of squares for
#' discrimination analysis.}
#'
#' @examples
#' Y <- sample(x = 0:1, size = 100, replace = TRUE, prob = NULL)
#' Ypred <- runif(n = 100)
#' result <- ConsensusOPLS:::DQ2(Ypred = Ypred, Y = Y)
#' result$DQ2
#' result$PRESSD
#' 
#' @keywords internal
#' @noRd
#' 
DQ2 <- function(Ypred, Y) {
    # Variable format control
    if (!all(Y %in% c(0, 1)))
        stop("Y must contain only 0 or 1 values.")
    if (length(Ypred) != length(Y))
        stop("Ypred and Y have different lengths.")

    # Residuals for each class
    E0 <- Ypred[Y==0]     #  - Y[Y==0]
    E1 <- Ypred[Y==1] - 1 #  - Y[Y==1]
    
    # Prediction error sum of squares for each class, disregarding the
    # prediction error when the class prediction is beyond the class label,
    # i.e. < 0 or > 1.
    PRESSD0 <- sum(E0[E0>0]^2)
    PRESSD1 <- sum(E1[E1<0]^2)
    
    # Total prediction error (PRESSD)
    PRESSD <- PRESSD0 + PRESSD1
    
    # Total Sum of Squares (TSS)
    TSS <- sum((Y-mean(Y))^2)

    # DQ2
    DQ2 <- 1 - (PRESSD / TSS)
    
    return (list("DQ2" = DQ2,
                 "PRESSD" = PRESSD))
}
