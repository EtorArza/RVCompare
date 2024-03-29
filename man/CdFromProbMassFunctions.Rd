% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RVCompare.R
\name{CdFromProbMassFunctions}
\alias{CdFromProbMassFunctions}
\title{The dominance rate of X_A over X_B for discrete distributions, given the probability mass functions.}
\usage{
CdFromProbMassFunctions(pMassA, pMassB)
}
\arguments{
\item{pMassA}{The probability mass function where pMassA[[i]] is the probability of x_i, p_A(x_i).}

\item{pMassB}{The probability mass function where pMassB[[i]] is the probability of x_i, p_B(x_i).}
}
\value{
Returns the dominance rate of X_A over X_B for discrete random variables.
}
\description{
Returns a real number in the interval [0,1] that represents the dominance rate of X_A over X_B.
}
\examples{
CdFromProbMassFunctions(c(0.2,0.6,0.2), c(0.3,0.3,0.4))
# > 0.6
# Notice how adding additional mass with the same cumulative distribution in both
# random variables does not change the result.
CdFromProbMassFunctions(c(0.2,0.6,0.2,0.2,0.2)/1.4, c(0.3,0.3,0.4,0.2,0.2)/1.4)
# > 0.6
}
