\docType{methods}
\name{designToJSON}
\alias{designToJSON}
\title{designToJSON}
\usage{
designToJSON(obj, expandX = FALSE)
}
\arguments{
  \item{design}{An object describing a study design}

  \item{expandX}{If true, the JSON representation will
  include a complete X matrix.}
}
\value{
The JSON representation of the design object
}
\description{
Generate a character string with the JSON representation of
a study design object, optionally with a complete X matrix.
Used to pass complex R objects into the supporting Java
library for the \code{rPowerlib} package.
}
\seealso{
\code{\link{design.glmmFG}} and \code{\link{design.glmmF}}
}

