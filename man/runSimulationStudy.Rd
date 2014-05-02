\name{runSimulationStudy}
\alias{runSimulationStudy}
\title{runSimulationStudy}
\usage{
runSimulationStudy(study.seed = 7634, study.data.dir = ".",
  study.figures.dir = ".", study.runEmpirical = TRUE)
}
\arguments{
  \item{study.seed}{the random number seed (defaults to
  7634)}

  \item{study.data.dir}{the directory into which data files
  are written (defaults to current working directory)}

  \item{study.figures.dir}{the directory into which pdf
  figures are written (defaults to the current working
  directory)}

  \item{study.runEmpirical}{if true, empirical power values
  will be recalculated. If false, existing empirical power
  values will be loaded from the R package.}
}
\description{
This function reproduces the simulation study results for
the manuscript:\cr Kreidler, S. M., Muller, K. E., &
Glueck, D. H. Calculating Power for the General Linear
Multivariate Model With One or More Gaussian Covariates, In
review.
}
\note{
The empirical power calculations may take several hours to
run
}

