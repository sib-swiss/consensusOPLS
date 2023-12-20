#' @title DQ2
#' @description
#' Calculates the discriminant Q2 value adjusted for values that are larger than
#' 1 for class 1 samples and smaller than 0 for 0 class samples. The idea is 
#' these values should not be penalized.
#' NB: It is assumed Y consists of ones and zeros indicating two classes.
#'
#' @param Ypred matrix. Predicted Y values.
#' @param Y matrix. Class labels (0 and 1).
#'
#' @return
#' \item{dqq}{ numeric. Discriminant Q2 value.}
#' \item{PRESSD}{ matrix. The total squared error for discrimination analysis.}
#'
#' @examples
#' Y <- sample(x = 0:1, size = 100, replace = TRUE, prob = NULL)
#' Ypred <- stats::runif(n = 100)  
#' result <-  ConsensusOPLS:::DQ2(Ypred = Ypred, Y = Y)
#' result$dqq
#' result$PRESSD
#' 
#' @keywords internal  

DQ2 <- function(Ypred, Y){
  # Variable format control
  if(!is.matrix(Ypred)){stop("Ypred is not a matrix.")}
  if(!is.matrix(Y)){
    stop("Y is not  matrix.")
  } else{
    test <- Y %in% c(0, 1)
    if(all(test) != TRUE){
      stop("Y must contain only 0 or 1 values.")
    }
  }
  
  # Find indices belonging to each Class
  Class0 <- which(Y == 0)
  Class1 <- which(Y == 1)
  
  # Calculate residuals for each Class
  E0 <- Ypred[Class0] - Y[Class0]   
  E1 <- Ypred[Class1] - Y[Class1]   
  
  # Find predictions for each Class
  E0count <- which(E0 > 0) # larger than 0
  E1count <- which(E1 < 0) # smaller than 0
  
  # Calculate SSE for each Class
  SSE0 <- sum(E0[E0count]**2)
  SSE1 <- sum(E1[E1count]**2)
  
  # Calculate total SSE (PRESSD)
  PRESSD <- SSE0 + SSE1             
  
  # Calculate Total Sum of Squares (TSS)
  Ym <- Y - mean(Y)                 
  TSS <- sum(Ym**2)
  
  # Return DQ2
  dqq <- 1 - (PRESSD / TSS)
  return(list("dqq" = dqq,
              "PRESSD" = PRESSD))
}
