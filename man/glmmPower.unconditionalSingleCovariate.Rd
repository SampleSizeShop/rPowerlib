\name{glmmPower.unconditionalSingleCovariate}
\alias{glmmPower.unconditionalSingleCovariate}
\title{glmmPower.unconditionalSingleCovariate}
\usage{
glmmPower.unconditionalSingleCovariate(design, hypothesis)
}
\arguments{
  \item{design}{a \code{design.glmmFG} object describing
  the study design}

  \item{hypothesis}{a \code{glh} object describing the
  hypothesis}
}
\value{
approximate power for the specified design and hypothesis
}
\description{
Calculate power for a multiple covariate design by
selecting the single most strongly correlated covariate and
applying the methods of:\cr Glueck, D. H., & Muller, K. E.
(2003). Adjusting power for a baseline covariate in linear
models. Statistics in Medicine, 22(16), 2535–2551.
}

