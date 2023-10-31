#' DQ2
#' Calculates the discriminant Q2 value adjusted for values that are larger than
#' 1 for class 1 samples and smaller than 0 for 0 class samples. The idea is 
#' these values should not be penalized.
#' NB: It is assumed Y consists of ones and zeros indicating two classes.
#'
#' @param Ypred: a matrix containing predicted Y values.
#' @param Y: a matrix containing class labels (0 and 1).
#'
#' @return
#' `dqq`: numeric. Discriminant Q2 value.
#' `PRESSD`: matrix. The total squared error for discrimination analysis.
#'
#' @examples
#' Y <- base::sample(0:1, 100, replace = TRUE)
#' Ypred <- stats::runif(100)  
#' result <- DQ2(Ypred, Y)
#' base::cat("DQ2:", result$dqq, "\n")
#' base::cat("PRESSD:", result$PRESSD, "\n")

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
  Class0 <- base::which(Y == 0)
  Class1 <- base::which(Y == 1)
  
  # Calculate residuals for each Class
  E0 <- Ypred[Class0] - Y[Class0]   
  E1 <- Ypred[Class1] - Y[Class1]   
  
  # Find predictions for each Class
  E0count <- base::which(E0 > 0) # larger than 0
  E1count <- base::which(E1 < 0) # smaller than 0
  
  # Calculate SSE for each Class
  SSE0 <- base::sum(E0[E0count]**2)
  SSE1 <- base::sum(E1[E1count]**2)
  
  # Calculate total SSE (PRESSD)
  PRESSD <- SSE0 + SSE1             
  
  # Calculate Total Sum of Squares (TSS)
  Ym <- Y - base::mean(Y)                 
  TSS <- base::sum(Ym**2)
  
  # Return DQ2
  dqq <- 1 - (PRESSD / TSS)
  return(list("dqq" = dqq,
              "PRESSD" = PRESSD))
}
