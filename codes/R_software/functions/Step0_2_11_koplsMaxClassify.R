#' @title koplsMaxClassify
#' @description Classification function that assesses class belonging of 
#' \code{data} based on the maximum value.
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
#' @param X matrix. It contains the predicted response matrix Y, where 
#' columns denote classes and rows observations.
#'
#' @return
#' \item{predClass}{ matrix. It contains the predicted class(es) of \code{data}.}
#'
#' @examples
#' data <- as.matrix(data.frame(x1 = c(2, 1, 5, 1),
#'                              x2 = c(7, 1, 1, 5),
#'                              x3 = c(9, 5, 4, 9),
#'                              x4 = c(3, 4, 1, 2)))
#' test <- koplsMaxClassify(X = data)
#' test
#' 
#' @keywords internal

koplsMaxClassify <- function(X){
  # Variable format control
  if(!is.matrix(X)){stop("X is not a matrix.")}
  
  # Search max position
  predClass <- base::matrix(data = base::apply(X = X, MARGIN = 1,
                                               FUN = function(row){
                                                 base::which.max(row)}),
                            ncol = 1)
  
  # Return the predicted class(es)
  return(predClass)
}
