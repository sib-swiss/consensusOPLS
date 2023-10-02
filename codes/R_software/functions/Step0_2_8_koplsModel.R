#' koplsModel
#' Function for training a K-OPLS model. The function constructs a
#' predictive regression model for predicting the values of 'Y' by
#' using the information in 'K'. The explained variation is separated
#' into predictive components (dimensionality is determined by the
#' parameter 'A') and 'Y'-orthogonal components (dimensionality determined 
#' by the parameter 'nox').
#'
#' @param K: Kernel matrix (un-centered); K = <phi(Xtr),phi(Xtr)>.
#' @param Y: Response matrix (un-centered/scaled). 
#' @param A: Number of predictive components. 
#' @param nox: Number of Y-orthogonal components. 
#' @param preProcK: Pre-processing parameters for the 'K' matrix:
#' `mc` for mean-centering, `no` for no centering.  
#' @param preProcY: Pre-processing parameters for the 'Y' matrix:
#' `mc` for mean-centering, `uv` for mc + scaling to unit variance,
#' `pa` for mc + Pareto, 'no' for no scaling. 
#'
#' @return
#'  `model` = a list with the following entries:
#'  `Cp` = Y loading matrix.
#'  `Sp` = Sigma matrix, containing singular values from Y'*K*Y used for scaling.
#'  `Up` = Y score matrix.
#'  `Tp` = Predictive score matrix for all Y-orthogonal components.
#'  `T` = Predictive score matrix for the final model.
#'  `co` = Y-orthogonal loading vectors.
#'  `so` = Eigenvalues from estimation of Y-orthogonal loading vectors.
#'  `To` = Y-orthogonal score matrix.
#'  `toNorm` = Norm of the Y-orthogonal score matrix prior to scaling.
#'  `Bt` = T-U regression coefficients for predictions.
#'  `A` = Number of predictive components.
#'  `nox` = Number of Y-orthogonal components.
#'  `K` = The kernel matrix.
#'  `EEprime` = The deflated kernel matrix for residual statistics.
#'  `sstot_K` = Total sums of squares in 'K'.
#'  `R2X` = Cumulative explained variation for all model components.
#'  `R2XO` = Cumulative explained variation for Y-orthogonal model components.
#'  `R2XC` = Explained variation for predictive model components after
#'  addition of Y-orthogonal model components.
#'  `sstot_Y` = Total sums of squares in Y.
#'  `R2Y` = Explained variation of Y.
#'  `preProcK` = Pre-processing setting for K.
#'  `preProcY` = Pre-processing setting for Y.
#'  `preProcParamsY` = Pre-processing scaling parameters for Y.
#'
#' @examples
#' K <- matrix(1:25, nrow = 5)
#' uniqueClass <- unique(TRUE)
#' preProcK <- "mc"
#' Y <- matrix(1:15, nrow = 5)
#' preProcY <- "mc"
#' A <- 2
#' nox <- 4
#' 
koplsModel <- function(K, Y, A, nox, preProcK, preProcY){
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
  
  # Initialize parameters
  I <- base::diag(length(K[, 1]))
  Kmc <- K
  
  # Check kernel centering
  if(preProcK == "mc"){
    Kmc <- koplsCenterKTrTr(K)
  }
  K <- base::vector("list", length = nox+1)
  K[[1]] <- Kmc
  
  # Save a copy of Y
  Y_old <- Y
  
  # Preprocess Y
  if(preProcY != "no"){
    if(preProcY == "mc"){
      scale <- "no"
    } else{
      scale <- preProcY
    }
    scaleParams <- koplsScale(Y, "mc", scale)
    Y <- scaleParams$matrix
  }
  
  # KOPLS model estimation
  ## step 1: SVD of Y'KY
  CSV <- base::svd(t(Y) %*% K[[1]] %*% Y)
  Cp <- CSV$u[, 1:A]
  Sp <- diag(CSV$d[1:A])
  
  ## step 2: Define Up
  Up <- Y %*% Cp
  
  ## step3: Loop over nox iterations
  Tp <- c() ; Bt <- c() ; co <- c()
  to <- c() ; toNorm <- c()
  for(i in 1:nox){
    ## step 4: Compute Tp
    Tp <- base::append(Tp, t(K[[1]][, i]) %*% Up %*% solve(Sp**(1/2)) )
    Bt <- base::append(Bt, solve( t(Tp[i])%*%Tp[i] ) %*% t(Tp[i]) %*% Up)
    
    ## step 5: SVD of T'KT
    temp <- base::svd( t(Tp[i])%*% (K[i,i]-Tp[i]%*%t(Tp[i])) %*% Tp[i]) 
    co <- base::append(co, temp$u[, 1])
    so <- base::append(so, temp$d[1])
    
    ## step 6: to
    to <- base::append(to, 
                       (K[i,i]-Tp[i]%*%t(Tp[i])) %*%Tp[i]%*%co[i]%*%(so[i]**{-1/2}))
    
    ## step 7: toNorm
    toNorm <- base::append(toNorm, sqrt(t(to[i])%*%to[i]) )
    
    ## step 8: Normalize to
    to[i] <- to[i]/toNorm[i]
    
    ## step 9: Update K
    K[1, i+1] <- K[1,i]%*%(I-to[i]%*%t(to[i]))
    
    ## step 10: Update Kii
    K[i+1, i+1] <- (I-to[i]%*%t(to[i]))%*%K[i, i]%*%(I-to[i]%*%t(to[i]))
  } ## step 11: end loop
  
  ## step 12: Tp[[nox+1]]
  Tp[nox+1] <- t(K[1, nox+1])%*%Up%*%(Sp**(-1/2))
  
  ## step 13: Bt[[nox+1]]
  Bt[i+1] <- solve( t(Tp[nox+1])%*%Tp[nox+1] )%*%t(Tp[nox+1])%*%Up
  
  # ---------- extra stuff -----------------
  # should work but not fully tested (MB 2007-02-19)
  sstot_Y = sum(sum(Y**2));
  F=Y-Up*t(Cp)
  R2Y = 1 - sum(sum( F**2 ))/sstot_Y
  # --------- #
  
  EEprime <- K[nox+1, nox+1]-Tp[nox+1]%*%t(Tp[nox+1])
  sstot_K <- sum(base::diag(K[1,1]))
  
  R2X <- c() ; R2X0 <- c() ; R2XC <- c() ; R2XO <- c() ; R2Yhat <- c()
  for(i in 1:(nox+1)){
    rss <- sum( base::diag(K[i,i]-Tp[i]%*%t(Tp[i])) )
    R2X <- base::append(R2X, 1- rss/sstot_K)
    
    rssc <- sum( base::diag(K[1,1]-Tp[i]%*%t(Tp[i])) )
    R2XC <- base::append(R2XC, 1- rssc/sstot_K)
    
    rsso = sum( base::diag( K[i,i] ))    
    R2XO = base::apprend(R2XO, 1 - rsso/sstot_K)
    
    # R2Yhat 22 Jan 2010 / MR - not fully tested
    Yhat <- Tp[i] %*% Bt[i] %*% t(Cp)
    R2Yhat <- base::append(R2Yhat, 1 - sum((Yhat - Y)**2)/sstot_Y )
  } # fin K-OPLS model
  
  # Convert to matrix structure
  To <- base::matrix(0, nrow = nrow(Tp[nox+1]), ncol = nox)
  for(i in 1:length(to)){
    To[, i]<- to[i]
  }
  
  return(model = list("Cp" = Cp, 
                      "Sp" = Sp, 
                      "Up" = Up,
                      "Tp" = Tp, 
                      "T" = Tp[nox+1], 
                      "co" = co,
                      "so" = so, 
                      "toNorm" = toNorm, 
                      "Bt" = Bt, 
                      "A" = A, 
                      "nox" = nox, 
                      "K" = K, 
                      "To" = To, 
                      "EEprime" =EEprime, 
                      "sstot_K" = sstot_K, 
                      "R2X" = R2X, 
                      "R2XO" = R2XO, 
                      "R2XC" = R2XC, 
                      "sstot_Y" = sstot_Y, 
                      "R2Y" = R2Y,
                      "R2Yhat" = R2Yhat, # R2Yhat 22 Jan 2010 / MR
                      "preProc" = list("K" = preProcK, "Y" = preProcY, 
                                       "paramsY" = scaleParams),
                      "class" = "kopls"))
}