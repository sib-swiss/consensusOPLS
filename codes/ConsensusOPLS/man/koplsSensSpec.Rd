% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/koplsSensSpec.R
\name{koplsSensSpec}
\alias{koplsSensSpec}
\title{koplsSensSpec}
\usage{
koplsSensSpec(trueClass, predClass, labelClass = NULL)
}
\arguments{
\item{trueClass}{matrix. Row vector of true class assignments (template).}

\item{predClass}{matrix. Matrix (or row vector) of class assignments to be 
compared.}

\item{labelClass}{A vector of all potential classes. Default, NULL, labels found
in \code{trueClass} and \code{codeClass} are considered.}
}
\value{
\item{data}{ A data frame containing, for each class (on lines), the number 
of True Positives \code{TP}, True Negatives \code{TN}, False Positives 
\code{FP} and False Negatives \code{FN}. It also contains (still in columns), 
the number of subjects \code{N}, the sensitivity \code{sens} and the 
specificity \code{spec}. The last two columns correspond to the mean 
sensitivity \code{meanSens} and mean specificity \code{meanSpec} (which are 
equals to the value for each class). The last line shows the total values for 
each column (named \code{tot}).}
}
\description{
Calculates sensitivity and specificity in a class-wise fashion.
}
\examples{
trueClass <- sample(x = 2:6, size = 100, replace = TRUE)
predClass <- sample(x = 1:5, size = 100, replace = TRUE)
test <- ConsensusOPLS:::koplsSensSpec(trueClass = trueClass, 
                                      predClass = predClass)
test

}
\keyword{internal}
