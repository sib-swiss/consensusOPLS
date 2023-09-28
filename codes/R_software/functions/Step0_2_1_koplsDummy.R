
#' KoplsDummy
#' Transform a vector into a dummy (binary) vector
#'
#' @param class: integer vector with classes to define.
#' @param numClasses: pre-defined the number of classes in the output. By 
#' default is equal to the length of `class` vector.
#'
#' @return
#' matrix: A matrix with rows corresponding to observations and columns 
#' to classes. Each element in matrix is either one (observation belongs to 
#' class) or zero (observation does not belong to class).
#' labels_sorted: The class labels that are found in class in sorted order.
#'
#' @examples
#' class <- c(5, 1, 2, 3, 4, 3, 2, 4, 3, 1, 3)
#' numClasses <- 2

koplsDummy <- function(class, numClasses = NA){
  # Define parameters
  labels <- base::unique(class)
  labels_sorted <- base::sort(labels, decreasing = FALSE)
  samples_number <- base::length(class)
  
  # Variable format control
  if(is.na(class)){stop("class must be contain an integer vector.")}
  if(is.na(numClasses)){
    sample_labels <- base::length(labels)
    ncol <- sample_labels
    default <- TRUE
  } else{
    sample_labels <- base::cut(x = class, breaks = numClasses,
                               ordered_result = TRUE)
    ncol <- base::length(levels(sample_labels))
    default <- FALSE
  }
  
  # Matrix initialization
  dummy <- base::matrix(data = 0, nrow = samples_number, ncol = ncol)
  
  # Search for class membership
  for(i in 1:ncol){
    if(default){
      dummy[class %in% labels_sorted[i], i] <- 1
    } else{
      labels_cut <- cbind(lower <- as.numeric(base::sub("\\((.+),.*", "\\1", 
                                                        levels(sample_labels)[i],
                                                        fixed = FALSE)),
                          upper <- as.numeric(base::sub("[^,]*,([^]]*)\\]", "\\1", 
                                                        levels(sample_labels)[i],
                                                        fixed = FALSE)))
      dummy[(class > labels_cut[1, "lower"] & 
               class <= labels_cut[1, "upper"]), i] <- 1
    }
  }
  # Return a list with 2 elements
  return(list("matrix" = dummy,
              "labels_sorted" = labels_sorted))
}