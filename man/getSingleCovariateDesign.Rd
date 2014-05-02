\name{getSingleCovariateDesign}
\alias{getSingleCovariateDesign}
\title{Create a design with a single covariate from a full GLMM(F,G) design}
\usage{
getSingleCovariateDesign(design)
}
\arguments{
  \item{design}{a design.glmmF or design.glmmFG object}
}
\value{
design.glmmFG object with a single covariate
}
\description{
Creates a new design.glmmFG object with a single covariate.
The covariate with the strongest covariance with any of the
outcomes is retained in the design
}

