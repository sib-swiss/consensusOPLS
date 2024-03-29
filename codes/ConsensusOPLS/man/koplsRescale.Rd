% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/koplsScale.R
\name{koplsRescale}
\alias{koplsRescale}
\title{koplsRescale}
\usage{
koplsRescale(scaleS, varargin = NULL)
}
\arguments{
\item{scaleS}{list. It contains scaling parameters.}

\item{varargin}{matrix. If defined, this matrix will be scaled and returned.
Otherwise the original data set in the \code{scaleS} object will be scaled 
and returned.}
}
\value{
A list containing the following entries:
\item{centerType}{ character. Indicates the centering method of the X matrix:
\code{mc} (mean-centering) or \code{no} (no centering).}
\item{scaleType}{ character. Indicates the scaling method of the X matrix: 
\code{uv} (unit variance), \code{pa} (pareto) or \code{no} (no scaling).}
\item{meanV}{ vector. Contains the mean values for all columns in X.}
\item{stdV}{ vector. Contains the standard deviations for all columns in X.}
\item{X}{ matrix. Scaled version of \code{varargin}, if defined. 
Otherwise, scaled version of \code{scaleS$X} from input. Scaling is done 
according to \code{centerType} and \code{scaleType}.
}
}
\description{
Scales a matrix based on pre-defined parameters from a scaling
object defined in a list (result of \code{koplsScale} function).
}
\examples{
data <- matrix(data = c(-1.732051, 0, 1.732051, 
                        2, 0,-2,
                        -1.732051, 0, 1.732051),
               nrow = 3, ncol = 3)
scaleS <- list("centerType" = "mc", "scaleType" = "pa", "meanV" = 0, 
               "stdV" = 1.581139, "matrix" = data)
test <- ConsensusOPLS:::koplsRescale(scaleS = scaleS, varargin = NULL)
test
test$X

}
\keyword{internal}
