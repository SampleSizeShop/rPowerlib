\docType{methods}
\name{simulateData}
\alias{simulateData}
\title{simulateData}
\usage{
simulateData(design, replicates = 10000, blockSize = 1000,
  outputDir = ".", filePrefix = "simulatedData", xNames = NA,
  yNames = NA, ...)
}
\arguments{
  \item{design}{An object describing a study design,
  commonly a \code{design.glmmFG}}

  \item{replicates}{the total number of data sets to
  generate}

  \item{blockSize}{the data sets can be written to multiple
  files.  The \code{blockSize} is the number of data sets
  which will be contained in each file.}

  \item{outputDir}{the directory to which the data set
  files will be written}

  \item{filePrefix}{the filename prefix for each data set
  file}

  \item{xNames}{the column names for the predictors in the
  data set}

  \item{yNames}{the column names for the outcomes in the
  data set}
}
\value{
Writes one or more CSV files containing data sets to the
output directory.
}
\description{
Simulate data sets based on the specified design.
}
\seealso{
\code{\link{design.glmmFG}} and \code{\link{design.glmmF}}
}

