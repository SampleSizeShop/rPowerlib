\docType{methods}
\name{simulateXMatrix}
\alias{simulateXMatrix}
\title{simulateXMatrix}
\usage{
simulateXMatrix(design)
}
\arguments{
  \item{design}{An object describing a study design,
  commonly a \code{design.glmmFG}}
}
\value{
A complete design matrix
}
\description{
Generate a complete X matrix for a given design.  If random
covariates are present in the design, the values of the
covariates will be randomly generated.
}
\seealso{
\code{\link{design.glmmFG}}
}

