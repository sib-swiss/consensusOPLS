% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/koplsPredict.R
\name{koplsPredict}
\alias{koplsPredict}
\title{koplsPredict}
\usage{
koplsPredict(KteTr, Ktest, Ktrain, model, nox, rescaleY = FALSE)
}
\arguments{
\item{KteTr}{matrix. The hybrid test/training kernel matrix. 
KteTr = <phi(Xte),phi(Xtr)>.}

\item{Ktest}{matrix. The pure test kernel matrix. Ktest = <phi(Xte),phi(Xte)>.}

\item{Ktrain}{matrix. The training kernel matrix (same as used in model 
training). Ktrain = <phi(Xtr),phi(Xtr)>.}

\item{model}{list. K-OPLS model object.}

\item{nox}{numeric. Number of Y-orthogonal components. If not specified, the 
number used during model training will be employed.}

\item{rescaleY}{logical. If \code{TRUE}, predicted values of the
response (Yhat) is rescaled according to the pre-processing settings of 
the model. If \code{FALSE}, Yhat is not rescaled (default).

#TODO: Why Tp, to, T, ... are produced here?}
}
\value{
A list with the following entries:
\item{Tp}{ matrix. Predicted predictive score matrix for all generations 0: 
\code{nox} of Y-orthogonal vectors.}
\item{to}{ vector. Predicted Y-orthogonal score vectors.}
\item{T}{ matrix. Predictive score matrix for the final model.}
\item{KteTr}{ matrix. Predictive score matrix for the final model with 
\code{nox} Y-orthogonal vectors.}
\item{EEprime}{ matrix. Calculated residuals for the test kernel \code{Ktest}, 
useful e.g. for residual statistics.}
\item{Yhat}{ matrix. Predicted values of the response matrix.}
}
\description{
Performs prediction of new samples from an existing K-OPLS model.
The function projects the Y-predictive and Y-orthogonal scores components 
to predict a value of the response matrix Y. The dimensions of the 
parameters is determined from the specified model.
}
\examples{
Xte <- matrix(data = stats::rnorm(n = 20), ncol=5)
Xtr <- matrix(data = stats::rnorm(n = 25), ncol=5)
KteTe <- ConsensusOPLS:::koplsKernel(X1 = Xte, X2 = Xte, 
                                     type='p', params=c(order=1.0))
KteTr <- ConsensusOPLS:::koplsKernel(X1 = Xte, X2 = Xtr, 
                                     type='p', params=c(order=1.0))
KtrTr <- ConsensusOPLS:::koplsKernel(X1 = Xtr, X2 = Xtr, 
                                     type='p', params=c(order=1.0))

Y <- matrix(data = stats::rnorm(n = 5), nrow = 5)
A <- 2
nox <- 4
preProcK <- "mc"
preProcY <- "mc"
model <- ConsensusOPLS:::koplsModel(K = KtrTr, Y = Y, A = A, nox = nox, 
                                    preProcK = preProcK, preProcY = preProcY)
pred <- ConsensusOPLS:::koplsPredict(KteTr = KteTr, Ktest = KteTe, 
                                     Ktrain = KtrTr, model = model, nox = nox,
                                     rescaleY = FALSE)
pred
}
\keyword{internal}
