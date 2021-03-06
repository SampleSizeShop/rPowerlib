% Generated by roxygen2 (4.0.2): do not edit by hand
\docType{class}
\name{glh}
\alias{glh}
\title{glh}
\description{
Class describing a general linear hypothesis for a multivariate model
}
\section{Slots}{

\describe{
\item{\code{alpha}}{The Type I error rate.  Must be a \code{numeric} value between 0 and 1.}

\item{\code{betweenContrast}}{The \code{matrix} of between participant contrasts.}

\item{\code{withinContrast}}{The \code{matrix} of within participant contrasts.}

\item{\code{thetaNull}}{The \code{matrix} of between participant contrasts.}

\item{\code{test}}{A \code{character} string indicating the statistical test. At present,
only the Hotelling-Lawley trace is supported with value "Hotelling-Lawley".}
}}
\note{
For theoretical details, please see

Muller, K. E., Lavange, L. M., Ramey, S. L., & Ramey, C. T. (1992).
Power Calculations for General Linear Multivariate Models Including
Repeated Measures Applications. Journal of the American Statistical
Association, 87(420), 1209-1226.
}

