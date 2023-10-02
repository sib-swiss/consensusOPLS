#' koplsMaxClassify
#' Classification function that assesses class belonging of 'data' based on 
#' the maximum value.
#'
#' @param data: matrix containing the predicted response matrix Y, where 
#' columns denote classes and rows observations.
#'
#' @return
#' `predClass`: the predicted class(es) of `data`.
#'
#' @examples
#'data <- as.matrix(data.frame(x1 = c(2, 1, 5, 1),
#'                             x2 = c(7, 1, 1, 5),
#'                             x3 = c(9, 5, 4, 9),
#'                             x4 = c(3, 4, 1, 2)))
#'test <- koplsMaxClassify(data = data)
#'test

koplsMaxClassify <- function(data){
  # Variable format control
  if(!is.matrix(data)){stop("data is not a matrix.")}
  
  # Initialization
  predClass <- numeric(nrow(data))
  
  # Search max position
  for(i in 1:nrow(data)){
    tmp <- which(data[i, ] == max(data[i, ]))
    if (length(tmp) == 1) {
      predClass[i] <- tmp
    } else {
      predClass[i] <- NaN
    }
  }
  
  # Transpose
  predClass <- t(as.vector(predClass))
  
  # Return the predicted class(es) of data
  return(predClass)
}
