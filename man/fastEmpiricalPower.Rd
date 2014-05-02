\docType{methods}
\name{fastEmpiricalPower}
\alias{fastEmpiricalPower}
\title{fastEmpiricalPower}
\usage{
fastEmpiricalPower(design, glh, realizations = 1000, replicates = 1000)
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
Calculate empiricalPower for a given design and hypothesis.
Utilizes Java code to provide faster computation of
empirical power.
}
\seealso{
\code{\link{design.glmmFG}} and \code{\link{glh}}
}

