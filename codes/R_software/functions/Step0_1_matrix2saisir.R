#' matrix2saisir
#' transforms a matrix in a `saisir` structure.
#' `Saisir means "statistique appliquée à l'interpretation des spectres 
#' infrarouge"
#' or "statistics applied to the interpretation of IR spectra".
#' 
#' @param data: R matrix
#' @param coderow (optional) : character corresponding to a code added 
#' before to the row identifiers. By default is NULL.
#' @param codecol (optional): character corresponding to a code added before
#' to the variables identifiers. By default is NULL.
#' @param position (optional): vector of character to indicate if the `coderow` 
#' and `codecol` is to be added `before`or `after` the number of the rows, or
#' respectively the variables identifiers. By default is `before`.
#'
#' @return data: SAISIR matrix which contains the data, the modified rownames,
#' and the modified colnames.

#' @examples
#' # create a matrix
#' A <- base::matrix(data = 1:8, nrow = 2, ncol = 4)
#' # apply the current function
#' B <- matrix2saisir(data = A, coderow = 'row #', codecol = 'Column #')
#' # explore all the elements contained in the results with:
#' head(B)
#' rownames(B)
#' colnames(B)

matrix2saisir <- function(data, 
                          coderow = NA, 
                          codecol = NA,
                          position_coderow = NA,
                          position_codecol = NA){
  # Default values
  if(is.na(coderow)){coderow <- ""}
  if(is.na(codecol)){codecol <- ""}
  if(is.na(position_coderow)){position_coderow <- c("before")}
  if(is.na(position_codecol)){position_codecol <- c("before")}
  
  # Variable format control 
  if(!is.numeric(data)){stop("data is not numeric.")}
  if(!is.character(coderow)){stop("coderow is not character.")}
  if(!is.character(codecol)){stop("codecol is not character.")}
  if(!is.character(position_coderow)){
    stop("position_coderow is not a character.")
  } else{
    if(position_coderow != "before" & position_coderow != "after"){
      stop("position_coderow must be `before` or `after`.")
    }
  }
  if(!is.character(position_codecol)){
    stop("position_codecol is not a character.")
  } else{
    if(position_codecol != "before" & position_codecol != "after"){
      stop("position_codecol must be `before` or `after`.")
    }
  }
  
  # Change the row and col names
  n <- nrow(data)
  p <- ncol(data)
  if(position_coderow == "before"){
    rownames(data) <- base::paste0(coderow, " ", base::as.character(1:n))
  } else{
    rownames(data) <- base::paste0(base::as.character(1:n), " ", coderow)
  }
  if(position_codecol == "before"){
    colnames(data) <- base::paste0(codecol, " ", base::as.character(1:p))
  } else{
    colnames(data) <- base::paste0(base::as.character(1:p), " ", codecol)
  }
  
  # Return the final matrix
  return(data)
}
