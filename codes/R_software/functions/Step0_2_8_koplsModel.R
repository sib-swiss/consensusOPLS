#' koplsModel
#' Function for training a K-OPLS model. The function constructs a
#' predictive regression model for predicting the values of 'Y' by
#' using the information in 'K'. The explained variation is separated
#' into predictive components (dimensionality is determined by the
#' parameter 'A') and 'Y'-orthogonal components (dimensionality determined 
#' by the parameter 'nox').
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
#' @param K: matrix. Kernel matrix (un-centered); K = <phi(Xtr),phi(Xtr)>.
#' @param Y: matrix. Response matrix (un-centered/scaled). 
#' @param A: numeric. Number of predictive components. Default is 1.
#' @param nox: numeric. Number of Y-orthogonal components. Default is 1.
#' @param preProcK: character. Pre-processing parameters for the 'K' matrix:
#' `mc` for mean-centering, `no` for no centering. Default is 'no'.
#' @param preProcY: character. Pre-processing parameters for the 'Y' matrix:
#' `mc` for mean-centering + no scaling, `uv` for mc + scaling to unit variance,
#' `pa` for mc + scaling to Pareto, `no` for no centering + no scaling. 
#' Default is 'no'.
#'
#' @return
#'  `model` is a list with the following entries:
#'  `Cp`: matrix. Y loading matrix.
#'  `Sp`: matrix. Sigma matrix, containing singular values from Y'*K*Y used 
#'  for scaling.
#'  `Sps`: matrix. Scaled Sigma matrix, containing scaled singular values.
#'  `Up`: matrix. Y score matrix.
#'  `Tp`: list. Predictive score matrix for all Y-orthogonal components.
#'  `T`: matrix. Predictive score matrix for the final model.
#'  `co`: list. Y-orthogonal loading vectors.
#'  `so`: list. Eigenvalues from estimation of Y-orthogonal loading vectors.
#'  `to`: list. Weight vector for the i-th latent component of the KOPLS model.
#'  `To`: matrix. Y-orthogonal score matrix.
#'  `toNorm`: list. Norm of the Y-orthogonal score matrix prior to scaling.
#'  `Bt`: list. T-U regression coefficients for predictions.
#'  `A`: numeric. Number of predictive components.
#'  `nox`: numeric. Number of Y-orthogonal components.
#'  `K`: matrix. The kernel matrix.
#'  `EEprime`: matrix. The deflated kernel matrix for residual statistics.
#'  `sstot_K`: numeric. Total sums of squares in 'K'.
#'  `R2X`: numeric. Cumulative explained variation for all model components.
#'  `R2XO`: numeric. Cumulative explained variation for Y-orthogonal model 
#'   components.
#'  `R2XC`: numeric. Explained variation for predictive model components after
#'   addition of Y-orthogonal model components.
#'  `sstot_Y`: numeric. Total sums of squares in Y.
#'  `R2Y`: numeric. Explained variation of Y.
#'  `R2Yhat`: numeric. Variance explained by the i-th latent component of the 
#'   model.
#'  `preProc$K`: character. Pre-processing setting for K.
#'  `preProc$Y`: character. Pre-processing setting for Y.
#'  `preProc$paramsY`: character or logical (if 'NA'). Pre-processing scaling 
#'   parameters for Y.
#'
#' @examples
#' K <- base::matrix(26:50, nrow = 5)
#' Y <- base::matrix(1:15, nrow = 5)
#' A <- 2
#' nox <- 4
#' preProcK <- "mc"
#' preProcY <- "mc"
#' test <- koplsModel(K = K, Y = Y, A = A, nox = nox, 
#'                    preProcK = preProcK, preProcY = preProcY)
#' ls(test)

koplsModel <- function(K, Y, A = 1, nox = 1, preProcK = "no", preProcY = "no"){
  # Variable format control
  if(!is.matrix(K)){stop("K is not a matrix.")}
  if(!is.matrix(Y)){stop("Y is not a matrix.")}
  if(!is.numeric(A)){stop("A is not a numeric.")}
  if(!is.numeric(nox)){stop("nox is not a numeric.")}
  if(!is.character(preProcK)){stop("preProcK is not a character.")
  } else{
    if(!(preProcK %in% c("mc", "no"))){stop("preProcK must be `mc` or `no`.")}
  }
  if(!is.character(preProcY)){stop("preProcY is not a character.")
  } else{ 
    if(!(preProcY %in% c("mc", "uv", "pa", "no"))){
      stop("preProcY must be `mc`, `uv`, `pa` or `no`.")}
  }
  
  # Function loading control
  if (!exists("koplsCenterKTrTr", mode = "function")) {
    warning("Remember to load the source code for the `koplsCenterKTrTr` function.")
  }
  if (!exists("koplsScale", mode = "function")) {
    warning("Remember to load the source code for the `koplsScale` function.")
  }
  
  # Initialize parameters
  I <- base::diag(ncol(K))
  
  # Check kernel centering
  if(preProcK == "mc"){
    Kmc <- koplsCenterKTrTr(K = K)
  } else{
    Kmc <- K
  }
  # K <- base::vector("list", length = nox+1)
  # K[[1]] <- Kmc
  K <- matrix(list(), ncol = nox+1, nrow = nox+1)
  K[1,1] <- list(Kmc)
  
  # Save a copy of Y
  Y_old <- Y
  
  # Preprocess Y
  if(preProcY != "no"){
    scaleParams <- koplsScale(X = Y_old, 
                              centerType = ifelse(preProcY == "mc", 
                                                  yes = "mc", no = "no"), 
                              scaleType = ifelse(preProcY == "mc", 
                                                 yes = "no", no = preProcY))
    Y <- scaleParams$matrix
  }
  
  # KOPLS model estimation
  ## step 1: SVD of Y'KY
  CSV <- base::svd(t(Y) %*% K[1,1][[1]] %*% Y,
                   nu = A, nv = A)
  # Extract left singular vectors
  Cp <- CSV$u
  # Extract the singular values
  if( A > 1){
    Sp <- base::diag(CSV$d[1:A])
    Sps <- base::diag(CSV$d[1:A]^(-1/2)) #scaled version
  } else{
    Sp <- CSV$d[1]
    Sps <- CSV$d[1]^(-1/2) #scaled version
  }
  
  ## step 2: Define Up
  Up <- Y %*% Cp
  
  # Initiate Yorth related variables
  to<-list(); co<-list(); so<-list(); toNorm<-list();
  Tp<-list(); Bt<-list();
  if(nox > 0){
    ## step3: Loop over nox iterations
    for(i in 1:nox){
      ## step 4: Compute Tp
      Tp[[i]] <- t(K[1,i][[1]]) %*% Up %*% Sps 
      Bt[[i]] <- base::solve( t(Tp[[i]])%*%Tp[[i]] ) %*% t(Tp[[i]]) %*% Up
      
      ## step 5: SVD of T'KT
      temp <- base::svd( t(Tp[[i]]) %*% (K[i,i][[1]]-Tp[[i]]%*%t(Tp[[i]])) %*% Tp[[i]],
                         nu = 1, nv = 1) 
      co[[i]] <- temp$u
      so[[i]] <- temp$d[1]
      
      ## step 6: to
      to[[i]] <- (K[i,i][[1]] - Tp[[i]]%*%t(Tp[[i]])) %*% Tp[[i]] %*% 
        co[[i]] %*% so[[i]]**(-1/2)
      
      ## step 7: toNorm
      toNorm[[i]] <- c(sqrt( t(to[[i]]) %*% to[[i]] ))
      
      ## step 8: Normalize to
      to[[i]] <- to[[i]] / toNorm[[i]]
      
      ## step 9: Update K
      scale_matrix <- I - to[[i]] %*% t(to[[i]])
      K[1, i+1][[1]] <- K[1,i][[1]] %*% scale_matrix
      
      ## step 10: Update Kii
      K[i+1, i+1][[1]] <- scale_matrix %*% K[i, i][[1]] %*% scale_matrix
      
    } ## step 11: end loop
  }
  
  ## step 12: Tp[[nox+1]]
  Tp[[nox+1]] <- t(K[1, nox+1][[1]]) %*% Up %*% Sps
  
  ## step 13: Bt[[nox+1]]
  Bt[[i+1]] <- base::solve( t(Tp[[nox+1]]) %*% Tp[[nox+1]] ) %*% t(Tp[[nox+1]]) %*% Up
  
  # ---------- extra stuff -----------------
  # should work but not fully tested (MB 2007-02-19)
  sstot_Y <- sum(sum(Y**2))
  F <- Y - Up %*% t(Cp)
  R2Y <- 1 - sum(sum( F**2 ))/sstot_Y
  # --------- #
  
  EEprime <- K[nox+1, nox+1][[1]] - Tp[[nox+1]] %*% t(Tp[[nox+1]])
  sstot_K <- sum(base::diag(K[1,1][[1]]))
  
  R2X <- c(); R2XO <- c(); R2XC <- c(); R2Yhat <- c();
  for(i in 1:(nox+1)){
    rss <- sum( base::diag(K[i,i][[1]] -  Tp[[i]]%*%t(Tp[[i]])) )
    R2X <- c(R2X, 1- rss/sstot_K)
    
    rssc <- sum( base::diag( K[1,1][[1]] - Tp[[i]]%*%t(Tp[[i]]) ) )
    R2XC <- c(R2XC, 1- rssc/sstot_K)
    
    rsso <- sum( base::diag( K[i,i][[1]] ))    
    R2XO <- c(R2XO, 1- rsso/sstot_K)
    
    # R2Yhat 22 Jan 2010 / MR - not fully tested
    Yhat <- Tp[[i]] %*% Bt[[i]] %*% t(Cp)
    R2Yhat <- c(R2Yhat, 1 - sum( sum((Yhat - Y)**2) )/sstot_Y )
  } # fin K-OPLS model
  
  # Convert to matrix structure
  if (nox > 0) {
    To <- base::matrix(base::unlist(to), nrow = nrow(Tp[[nox+1]]), ncol = nox, byrow=FALSE)
  } else {
    To <- NULL
  }
  
  return(list("Cp" = Cp, 
              "Sp" = Sp, 
              "Sps" = Sps, 
              "Up" = Up,
              "Tp" = Tp, 
              "T" = as.matrix(Tp[[nox+1]]), 
              "co" = co,
              "so" = so, 
              "to" = to, 
              "To" = To,
              "toNorm" = toNorm, 
              "Bt" = Bt, 
              "A" = A, 
              "nox" = nox, 
              "K" = K, 
              
              #extra stuff
              "EEprime" =EEprime, 
              "sstot_K" = sstot_K, 
              "R2X" = R2X, 
              "R2XO" = R2XO, 
              "R2XC" = R2XC, 
              "sstot_Y" = sstot_Y, 
              "R2Y" = R2Y,
              "R2Yhat" = R2Yhat, # R2Yhat 22 Jan 2010 / MR
              
              #pre-processing
              "preProc" = list("K" = preProcK, 
                               "Y" = preProcY, 
                               "paramsY" = ifelse(test = (preProcY != "no"),
                                                  yes = scaleParams,
                                                  no = NA)),
              "class" = "kopls"))
}
