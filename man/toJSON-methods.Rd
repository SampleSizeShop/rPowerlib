\docType{methods}
\name{toJSON}
\alias{toJSON}
\alias{toJSON,glh,glh-method}
\alias{toJSON,glh-method}
\title{toJSON}
\usage{
toJSON(obj)

\S4method{toJSON}{glh}(obj)
}
\arguments{
  \item{design}{An object describing a study design}

  \item{expandX}{If true, the JSON representation will
  include a complete X matrix.}
}
\value{
The JSON representation of the object
}
\description{
Generate a character string with the JSON representation of
an object, most commonly a \code{glh} object. Used to pass
complex R objects into the supporting Java library for the
\code{rPowerlib} package.
}
\seealso{
\code{\link{glh}}
}

