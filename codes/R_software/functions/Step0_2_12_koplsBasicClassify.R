#' koplsBasicClassify
#' Classification function that assesses class belonging of a predicted
#' response in 'data' based on a fixed threshold 'k'.
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
#' @param data: matrix containing the predicted response matrix Y,
#' where columns denote classes and rows observations.
#' @param k: threshold value used to assign class categories.
#'
#' @return
#' `predClass`: the predicted class(es) of `data` given `k`.
#'
#' @examples
#' data <- as.matrix(data.frame(x1 = c(2, 1, 5, 1),
#'                              x2 = c(7, 1, 1, 5),
#'                              x3 = c(9, 5, 4, 9),
#'                              x4 = c(3, 4, 1, 2)))
#' k <- 4
#' test <- koplsBasicClassify(data = data, k = k)
#' test

koplsBasicClassify <- function(data, k){
  # Variable format control
  if(!is.matrix(data)){stop("data is not a matrix.")}
  if(!is.numeric(k)){stop("k is not numeric.")}
  
  # Initialization
  predClass <- base::matrix(NA, nrow = nrow(data), ncol = ncol(data))
  
  # Search predicted class(es)
  for (i in 1:nrow(data)) {
    tmp <- base::which(data[i, ] > k)
    predClass[i, 1:length(tmp)] <- tmp
  }
  
  # Return the predicted class(es) of 'data' given 'k'.
  return(predClass)
}
