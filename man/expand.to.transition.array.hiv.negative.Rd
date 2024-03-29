% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/transitions.R
\name{expand.to.transition.array.hiv.negative}
\alias{expand.to.transition.array.hiv.negative}
\title{Expands an array into a transition array for the HIV-negative population for a given dimension}
\usage{
expand.to.transition.array.hiv.negative(jheem, to.expand, transition.dimension)
}
\arguments{
\item{jheem}{An object of \code{\link{jheem-class}}}

\item{to.expand}{The array to expand. At minimum, must have .from and .to dimensions for the the transition dimension. May also have age, race, subpopulation, sex, risk, non.hiv.subset}

\item{transition.dimension}{The name of the dimension for which this is a transition array. Must be one of age, race, subpopulation, sex, risk, non.hiv.subset}
}
\description{
Expands an array into a transition array for the HIV-negative population for a given dimension
}
\seealso{
Other Managing transition arrays: 
\code{\link{expand.to.transition.array.hiv.positive}()}
}
\concept{Managing transition arrays}
