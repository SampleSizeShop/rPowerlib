\name{getFixedDesign}
\alias{getFixedDesign}
\title{Extract the fixed-only portion of a GLMM(F,G) design.}
\usage{
getFixedDesign(design)
}
\arguments{
  \item{design}{a design.glmmF or design.glmmFG object}
}
\value{
design.glmmF object
}
\description{
Creates a design.glmmF object from a design.glmmFG or
design.glmmF object, after removing information regarding
covariates.
}

