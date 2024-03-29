% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/koplsCenterK.R
\name{koplsCenterKTeTr}
\alias{koplsCenterKTeTr}
\title{koplsCenterKTeTr}
\usage{
koplsCenterKTeTr(KteTr, KtrTr)
}
\arguments{
\item{KteTr}{matrix. Contains the hybrid test/training kernel matrix, 
KteTr = <phi(Xte), phi(Xtr)>.}

\item{KtrTr}{matrix. Contains the training kernel matrix; 
Ktrain = <phi(Xtr), phi(Xtr)>.}
}
\value{
The centered kernel matrix.
}
\description{
Centering function for the hybrid test/training kernel, which
is constructed from the test matrix Xte and the training matrix
Xtr as KteTr = <phi(Xte), phi(Xtr)>. Requires additional
(un-centered) training kernel to estimate mean values.
}
\examples{
Xte <- matrix(data = stats::rnorm(n = 20), ncol=5)
Xtr <- matrix(data = stats::rnorm(n = 25), ncol=5)
KteTr <- ConsensusOPLS:::koplsKernel(X1 = Xte, X2 = Xtr, 
                                     type='p', params=c(order=1.0))
KtrTr <- ConsensusOPLS:::koplsKernel(X1 = Xtr, X2 = Xtr, 
                                     type='p', params=c(order=1.0))
ConsensusOPLS:::koplsCenterKTeTr(KteTr = KteTr, KtrTr = KtrTr)

}
\keyword{internal}
