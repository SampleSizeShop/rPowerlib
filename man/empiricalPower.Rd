\docType{methods}
\name{empiricalPower}
\alias{empiricalPower}
\title{empiricalPower}
\usage{
empiricalPower(design, glh, replicates = 1000, realizations = 1000)
}
\arguments{
  \item{design}{An object describing a study design}

  \item{glh}{An object describing the general linear
  hypothesis to be tested}
}
\value{
empirical unconditional power
}
\description{
Calculate empirical unconditional power for a given design
and hypothesis.  This function is implemented purely in R
and may run slowly for complex designs.
}
\seealso{
\code{\link{design.glmmFG}} and \code{\link{glh}}
}

