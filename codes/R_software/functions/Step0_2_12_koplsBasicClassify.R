#' @title koplsBasicClassify
#' @description Classification function that assesses class belonging of a 
#' predicted response in \code{data} based on a fixed threshold \code{k}.
#' 
#' # ------------------------------------------------------------------------ #
#' This file is part of the K-OPLS package, developed by Max Bylesjo, 
#' University of Umea, Judy Fonville and Mattias Rantalainen, Imperial College.
#' 
#' Copyright (c) 2007-2010 Max Bylesjo, Judy Fonville and Mattias Rantalainen 
#' 
#' This code has been extended and adapted under the terms of the GNU General 
#' Public License version 2 as published by the Free Software Foundation.
#' # ------------------------------------------------------------------------ #
#'
#' @param X matrix. It contains the predicted response matrix Y,
#' where columns denote classes and rows observations.
#' @param k numeric. Threshold value used to assign class categories.
#'
#' @return
#' \item{predClass}{ matrix. The predicted class(es) of \code{data} given 
#' \code{k}.}
#'
#' @examples
#' data <- as.matrix(data.frame(x1 = c(2, 1, 5, 1),
#'                              x2 = c(7, 1, 1, 5),
#'                              x3 = c(9, 5, 4, 9),
#'                              x4 = c(3, 4, 1, 2)))
#' k <- 4
#' test <- koplsBasicClassify(X = data, k = k)
#' test
#' 
#' @keywords internal

koplsBasicClassify <- function(X, k){
  # Variable format control
  if(!is.matrix(X)){stop("X is not a matrix.")}
  if(!is.numeric(k)){stop("k is not numeric.")}
  
  # Search predicted class(es)
  predClass <- base::apply(X > k, MARGIN = 1, 
                           FUN = function(row) {base::which(row)})
  
  # Return the predicted class(es).
  return(predClass)
}
