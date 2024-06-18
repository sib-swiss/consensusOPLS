#' @title DQ2
#' @description
#' Calculates the discriminant Q2 value adjusted for values that are greater than
#' 1 for 1-class samples and smaller than 0 for 0-class samples. The idea is 
#' these values should not be penalized.
## TOCHECK: Don't understand !!!!!!!
#'
#' @param Ypred A vector of predicted values.
#' @param Y A vector of two-class 0/1 labels.
#'
#' @returns
#' \item{dqq}{ numeric. Discriminant Q2 value.}
#' \item{PRESSD}{matrix. The total squared error for discrimination analysis.}
#'
#' @examples
#' Y <- sample(x = 0:1, size = 100, replace = TRUE, prob = NULL)
#' Ypred <- stats::runif(n = 100)
#' result <- ConsensusOPLS:::DQ2(Ypred = Ypred, Y = Y)
#' result$dqq
#' result$PRESSD
#' 
#' @keywords internal
#' @noRd
#' 
DQ2 <- function(Ypred, Y) {
    # Variable format control
    if (!all(Y %in% c(0, 1))) stop("Y must contain only 0 or 1 values.")

    # Find indices belonging to each Class
    class0 <- which(Y == 0)
    class1 <- which(Y == 1)
    
    # Calculate residuals for each Class
    E0 <- Ypred[class0] - Y[class0]
    E1 <- Ypred[class1] - Y[class1]
    
    # Find predictions for each Class
    E0count <- which(E0 > 0) # larger than 0
    E1count <- which(E1 < 0) # smaller than 0
    
    # Calculate SSE for each Class
    SSE0 <- sum(E0[E0count]^2)
    SSE1 <- sum(E1[E1count]^2)
    
    # Calculate total SSE (PRESSD)
    PRESSD <- SSE0 + SSE1
    
    # Calculate Total Sum of Squares (TSS)
    Ym <- Y - mean(Y)
    TSS <- sum(Ym^2)

    # DQ2
    dqq <- 1 - (PRESSD / TSS)
    
    return (list("dqq" = dqq,
                 "PRESSD" = PRESSD))
}
