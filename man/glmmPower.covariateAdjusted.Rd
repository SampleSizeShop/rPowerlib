\name{glmmPower.covariateAdjusted}
\alias{glmmPower.covariateAdjusted}
\title{glmmPower.covariateAdjusted}
\usage{
glmmPower.covariateAdjusted(design, hypothesis,
  mAdjust = mAdjust.expectedProjection)
}
\arguments{
  \item{design}{a \code{design.glmmFG} object describing
  the study design}

  \item{hypothesis}{a \code{glh} object describing the
  hypothesis}

  \item{mAdjust}{(optional) a function calculating the
  adjustment factor for the M matrix}
}
\value{
approximate power for the specified design and hypothesis
}
\description{
Calculate power using the covariate adjustment approach
proposed in the manuscript:\cr Kreidler, S. M., Muller, K.
E., & Glueck, D. H. Calculating Power for the General
Linear Multivariate Model With One or More Gaussian
Covariates, In review.
}

