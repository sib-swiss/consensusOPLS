#' @title koplsModel
#' @description Function for training a K-OPLS model. The function constructs a
#' predictive regression model for predicting the values of 'Y' by
#' using the information in 'K'. The explained variation is separated
#' into predictive components (dimensionality is determined by the
#' parameter 'A') and 'Y'-orthogonal components (dimensionality determined 
#' by the parameter 'nox').
#'
#' @param K matrix. Kernel matrix (un-centered); K = <phi(Xtr),phi(Xtr)>.
#' @param Y matrix. Response matrix (un-centered/scaled). 
#' @param A numeric. Number of predictive components. Default is 1.
#' @param nox numeric. Number of Y-orthogonal components. Default is 1.
#' @param preProcK character. Pre-processing parameters for the 'K' matrix:
#' \code{mc} for mean-centering, \code{no} for no centering. Default is 'no'.
#' @param preProcY: character. Pre-processing parameters for the 'Y' matrix:
#' \code{mc} for mean-centering + no scaling, \code{uv} for mc + scaling to unit 
#' variance, \code{pa} for mc + scaling to Pareto, \code{no} for no centering + 
#' no scaling. Default is \code{no}.
#'
#' @return A list with the following entries:
#' \item{Cp}{ matrix. Y loading matrix.}
#' \item{Sp}{ matrix. Sigma matrix, containing singular values from Y'*K*Y used 
#'  for scaling.}
#' \item{Sps}{ matrix. Scaled Sigma matrix, containing scaled singular values.}
#' \item{Up}{ matrix. Y score matrix.}
#' \item{Tp}{ list. Predictive score matrix for all Y-orthogonal components.}
#' \item{T}{ matrix. Predictive score matrix for the final model.}
#' \item{co}{ list. Y-orthogonal loading vectors.}
#' \item{so}{ list. Eigenvalues from estimation of Y-orthogonal loading vectors.}
#' \item{to}{ list. Weight vector for the i-th latent component of the KOPLS 
#' model.}
#' \item{To}{ matrix. Y-orthogonal score matrix.}
#' \item{toNorm}{ list. Norm of the Y-orthogonal score matrix prior to scaling.}
#' \item{Bt}{ list. T-U regression coefficients for predictions.}
#' \item{A}{ numeric. Number of predictive components.}
#' \item{nox}{ numeric. Number of Y-orthogonal components.}
#' \item{K}{ matrix. The kernel matrix.}
#' \item{EEprime}{ matrix. The deflated kernel matrix for residual statistics.}
#' \item{sstot_K}{ numeric. Total sums of squares in 'K'.}
#' \item{R2X}{ numeric. Cumulative explained variation for all model components.}
#' \item{R2XO}{ numeric. Cumulative explained variation for Y-orthogonal model 
#' components.}
#' \item{R2XC}{ numeric. Explained variation for predictive model components 
#' after addition of Y-orthogonal model components.}
#' \item{sstot_Y}{ numeric. Total sums of squares in Y.}
#' \item{R2Y}{ numeric. Explained variation of Y.}
#' \item{R2Yhat}{ numeric. Variance explained by the i-th latent component of 
#' the model.}
#' \item{preProc$K}{ character. Pre-processing setting for K.}
#' \item{preProc$Y}{ character. Pre-processing setting for Y.}
#' \item{preProc$paramsY}{ character. Pre-processing scaling parameters for Y.}
#'
#' @examples
#' K <- matrix(26:50, nrow = 5)
#' Y <- matrix(1:15, nrow = 5)
#' A <- 2
#' nox <- 4
#' preProcK <- "mc"
#' preProcY <- "mc"
#' test <- ConsensusOPLS:::koplsModel(K = K, Y = Y, A = A, nox = nox, 
#'                                    preProcK = preProcK, preProcY = preProcY)
#' ls(test)
#' 
#' @keywords internal

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
  
  # Initialize parameters
  I <- diag(ncol(K))
  
  # Check kernel centering
  if(preProcK == "mc"){
    Kmc <- ConsensusOPLS:::koplsCenterKTrTr(K = K)
  } else{
    Kmc <- K
  }
  K <- matrix(data = list(), ncol = nox+1, nrow = nox+1)
  K[1,1] <- list(Kmc)
  
  # Save a copy of Y
  Y_old <- Y
  
  # Preprocess Y
  if(preProcY != "no"){
    scaleParams <- ConsensusOPLS:::koplsScale(X = Y_old, 
                                              centerType = 
                                                ifelse(preProcY == "mc", 
                                                       yes = "mc", no = "no"), 
                                              scaleType = 
                                                ifelse(preProcY == "mc", 
                                                       yes = "no", no = preProcY))
    Y <- scaleParams$matrix
  }
  
  # KOPLS model estimation
  ## step 1: SVD of Y'KY
  CSV <- svd(x = crossprod(Y, crossprod(t(K[1,1][[1]]), Y)),
             nu = A, nv = A)
  # Extract left singular vectors
  Cp <- CSV$u
  # Extract the singular values
  if( A > 1){
    Sp <- diag(CSV$d[1:A])
    Sps <- diag(CSV$d[1:A]^(-1/2)) #scaled version
  } else{
    Sp <- CSV$d[1]
    Sps <- CSV$d[1]^(-1/2) #scaled version
  }
  
  ## step 2: Define Up
  Up <- crossprod(t(Y), Cp)
  
  # Initiate Yorth related variables
  to <- list(); co <- list(); so <- list(); toNorm <- list();
  Tp <- list(); Bt <- list();
  i <- 1
  
  ## step3: Loop over nox iterations
  while(i <= nox){
    ## step 4: Compute Tp
    Tp[[i]] <- crossprod(K[1,i][[1]], tcrossprod(Up, t(Sps)) )
    Bt[[i]] <- crossprod( t(solve( crossprod(Tp[[i]]) )), 
                          crossprod(Tp[[i]], Up))
    
    ## step 5: SVD of T'KT
    temp <- svd(x = 
                  crossprod(Tp[[i]], 
                            tcrossprod(K[i,i][[1]] -
                                         tcrossprod(Tp[[i]]), 
                                       t(Tp[[i]]))),
                nu = 1, nv = 1) 
    co[[i]] <- temp$u
    so[[i]] <- temp$d[1]
    
    ## step 6: to
    to[[i]] <- tcrossprod(tcrossprod(tcrossprod(K[i,i][[1]] - 
                                                  tcrossprod(Tp[[i]]), 
                                                t(Tp[[i]])),
                                     t(co[[i]])), t(so[[i]]**(-1/2)))
    
    ## step 7: toNorm
    toNorm[[i]] <- c(sqrt( crossprod(to[[i]]) ))
    
    ## step 8: Normalize to
    to[[i]] <- to[[i]] / toNorm[[i]]
    
    ## step 9: Update K
    scale_matrix <- I - tcrossprod(to[[i]])
    K[1, i+1][[1]] <- tcrossprod(K[1,i][[1]], t(scale_matrix))
    
    ## step 10: Update Kii
    K[i+1, i+1][[1]] <- tcrossprod(scale_matrix, tcrossprod(t(scale_matrix),
                                                            t(K[i, i][[1]])))
    
    # Update i
    i <- i + 1
  }## step 11: end loop
  
  ## step 12: Tp[[nox+1]]
  Tp[[nox+1]] <- crossprod(K[1, nox+1][[1]], crossprod(t(Up), Sps))
  
  ## step 13: Bt[[nox+1]]
  Bt[[nox+1]] <- crossprod( t(solve( crossprod(Tp[[nox+1]]) )),
                            crossprod(Tp[[nox+1]], Up))
  
  # ---------- extra stuff -----------------
  # should work but not fully tested (MB 2007-02-19)
  sstot_Y <- sum( sum(Y**2))
  F <- Y - tcrossprod(Up, Cp)
  R2Y <- 1 - sum( sum( F**2 ))/sstot_Y
  # --------- #
  
  EEprime <- K[nox+1, nox+1][[1]] - tcrossprod(Tp[[nox+1]])
  sstot_K <- sum( diag(K[1,1][[1]]))
  
  R2X <- c(); R2XO <- c(); R2XC <- c(); R2Yhat <- c();
  for(i in 1:(nox+1)){
    rss <- sum( diag(K[i,i][[1]] -  tcrossprod(Tp[[nox+1]])) )
    R2X <- c(R2X, 1- rss/sstot_K)
    
    rssc <- sum( diag( K[1,1][[1]] - tcrossprod(Tp[[nox+1]]) ) )
    R2XC <- c(R2XC, 1- rssc/sstot_K)
    
    rsso <- sum( diag( K[i,i][[1]] ))    
    R2XO <- c(R2XO, 1- rsso/sstot_K)
    
    # R2Yhat 22 Jan 2010 / MR - not fully tested
    Yhat <- crossprod(t(Tp[[i]]), tcrossprod(Bt[[i]], Cp))
    R2Yhat <- c(R2Yhat, 1 - sum( sum((Yhat - Y)**2) )/sstot_Y )
  } # fin K-OPLS model
  
  # Convert to matrix structure
  if (nox > 0) {
    To <- matrix(data = unlist(to), nrow = nrow(Tp[[nox+1]]), 
                 ncol = nox, byrow = FALSE)
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
                                                  no = "no")),
              "class" = "kopls"))
}