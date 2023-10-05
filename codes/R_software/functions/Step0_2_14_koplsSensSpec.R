#' koplsSensSpec
#' Calculates sensitivity and specificity in a class-wise fashion.
#'
#' @param v: row vector of true class assignments (template). 
#' @param m: matrix (or row vector) of class assignments to be compared.
#'
#' @return
#' a list containing:
#' `sensvec`: sensitivity for each class.
#' `specvec`: specificity for each class.
#' `classvec`: the class identifier corresponding to each column in `sensvec` 
#' and `specvec.`
#' `tot_sens`: total sensitivity.
#' `meanSens`: mean sensitivity over all classes.
#' `meanSpec`: mean specificity over all classes.
#'
#' @examples
#'

koplsSensSpec <- function(v, m){
  # Variable format control
  if(!is.vector(v)){stop("v is not a vector.")}
  if(!is.matrix(m) & !is.vector(m)){stop("m is neither a matrix nor a vector.")}
  
  # Check dimensions
  if(is.vector(v) & is.vector(m)){
    if(length(m) != length(v)){
      stop("Template vector (v) differs in length from the vector (m) to be compared.")
    }
    
    # Forces matrix conversion for the code below
    m <- as.matrix(m)
  } else{
    if(is.vector(v) & is.matrix(m)){
      if(nrow(m) != length(v)){
        stop("Template vector (v) differs in length from the matrix (m) to be compared.")
      }
    }
  }
  
  # Define classes to compare
  classes <- sort(unique(v))
  
  # 
  if(ncol(m) > 1){
    mtemp <- c()
    vtemp <- c()
    
    for(j in 1:nrow(m)){
      for(i in 1:ncol(m)){
        if(i == 1){ # 1st column always added
          vtemp <- c(vtemp, v[j])
          
          if(is.na(m[j,i])){
            mtemp <- c(mtemp, max(classes)+1)
          } else{mtemp <- c(mtemp, m[j,i])}
        } else{
          if(!is.na(m[j,i])){
            vtemp <- c(vtemp, v[j])
            mtemp <- c(mtemp, m[j,i])
          }
        }
      }
    }
    v <- vtemp
    m <- mtemp
  }
  
  row_indices <- 1:length(v)
  
  # Initialisation parameters
  ind <- which(!is.na(m))
  sensvec <- numeric(length(classes))
  specvec <- numeric(length(classes))
  classvec <- numeric(length(classes))
  
  for (i in 1:length(classes)) {
    class <- classes[i]
    rows_inclass <- row_indices[v[ind] == class]
    rows_notinclass <- row_indices[v[ind] != class]
    
    tp_rate <- sum(m[ind][rows_inclass] == class) / length(rows_inclass)
    tn_rate <- sum(m[ind][rows_notinclass] != class) / length(rows_notinclass)
    
    sensvec[i] <- tp_rate
    specvec[i] <- tn_rate
    classvec[i] <- class
  }
  
  #Return the list
  return(list("sensvec" = sensvec, 
              "specvec" = specvec, 
              "classvec" = classvec, 
              "tot_sens" = sum(v == m)/length(v),
              "meanSens" = mean(sensvec),
              "meanSpec" = mean(specvec)))
}