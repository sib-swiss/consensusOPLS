#' @title MBVIP
#' @description Calculate the VIP (Variable Importance in Projection) for each variable of a 
#' Consensus-OPLS model.
#' 
#' @param collection list. The collection list containing each block of data.
#' @param Y matrix. The response matrix to predict.
#' @param model list. ConsensusOPLS model.
#'
#' @return 
#' A table with the results:
#' `VIP` = sqrt(p*q/s), with:
#' `p` is the number of variables in each block
#' `q` is the explained variance of Y associated to each variable
#' `s` is the total Y variance explained by the model
#'
#' @examples
#' TO DO

MBVIP <- function(collection, Y, model){
  # Variable format control
  if(!is.list(collection)){stop("data is not a list.")}
  if(!is.matrix(Y)){stop("Y is not a matrix.")}
  if(!is.list(model)){stop("model is not a list.")}

  ntable <- nrow(model$model$Model$loadings)
  VIP <- list()
  
  for (table in 1:ntable) {
    # Dimensions of the data in the collection
    p <- ncol(collection[[table]])
    m <- nrow(model$model$Model$T)
    h <- ncol(model$model$Model$T)
    
    W <- list()
    Q <- list()
    
    for (i in 1:model$model$Model$A) {
      q <- (t(Y) %*% model$model$Model$T[, i]) %*% 
        ( t(model$model$Model$T[, i]) %*% model$model$Model$T[, i] )**(-1)
      u <- (Y %*% q) %*% (t(q) %*% q)**(-1)
      w <- (t(collection[[table]]) %*% u) %*% (t(u) %*% u)**(-1)
      w <- w / base::norm(w, type = "F")
      W[[i]] <- w
      Q[[i]] <- q
    }
    
    Q <- t(Q[[1:model$model$Model$A]])
    s <- base::diag( t(model$model$Model$T) %*% model$model$Model$T %*% Q %*% t(Q) )
    
    VIPtemp <- numeric(p)
    
    for (i in 1:p) {
      weight <- sapply(1:h, function(j) (W[[j]][i, j] / base::norm(x = W[[j]],
                                                                  type = "F"))^2)
      q <- t(s) %*% weight
      VIPtemp[i] <- sqrt(p * q / sum(s))
    }
    
    VIP[[table]] <- VIPtemp
  }

  # Return the result
  return(VIP)
}
