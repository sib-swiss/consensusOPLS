% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/VariableImportanceProjection.R
\name{VIP}
\alias{VIP}
\title{Variable Importance in Projection}
\usage{
VIP(data, Y, model = NULL, mc.cores = 1, ...)
}
\arguments{
\item{data}{A list of data blocks.}

\item{Y}{A vector, factor, dummy matrix or numerical matrix for the response.}

\item{model}{A ConsensusOPLS model Default, NULL, a model will be constructed.}

\item{mc.cores}{Number of cores for parallel computing. Default: 1.}

\item{...}{arguments to pass to \code{RVConsensusOPLS}}
}
\value{
A table with the results:
\code{VIP = sqrt(p*q/s)}, where
\code{p} is the number of variables in each block,
\code{q} the explained variance of Y associated to each variable, and 
\code{s} the total Y variance explained by the model.
}
\description{
Calculate the VIP (Variable Importance in Projection) for each variable 
in a ConsensusOPLS model.
}
\examples{
vip <- VIP(data=demo_3_Omics[c("MetaboData", "MicroData", "ProteoData")], 
           Y=demo_3_Omics$Y)
str(vip)
}
