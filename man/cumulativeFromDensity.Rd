% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RVCompare.R
\name{cumulativeFromDensity}
\alias{cumulativeFromDensity}
\title{Get the cumulative distribution function given the distribution function.}
\usage{
cumulativeFromDensity(densityX, xlims, sanityChecks = TRUE)
}
\arguments{
\item{densityX}{The probability density function.}

\item{xlims}{the domain of definition of the density function.}

\item{sanityChecks}{(optional parameter, default = TRUE) boolean value indicating wether to check if the density function is correctly defined.}
}
\value{
a callable function representing the cumulative distribution.
}
\description{
Get the cumulative distribution function given the distribution function.
}
\examples{
cumulativeProbability <- cumulativeFromDensity(normalDensity(0,1), c(-4,4), FALSE)
x <- seq(-4, 4, length.out=101)
plot(x, cumulativeProbability(x), type = "l")
}
\keyword{internal}
